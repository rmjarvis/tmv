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


#ifndef TMV_SmallMatrixArith_H
#define TMV_SmallMatrixArith_H

#include "TMV_SmallVectorArith.h"
#include "TMV_SmallMatrixArithFunc.h"

#define CT std::complex<T>
#define CCT ConjRef<std::complex<T> >
#define VCT VarConjRef<std::complex<T> >

namespace tmv {

  template <class T, size_t M, size_t N> class SmallMatrixComposite : 
    public GenSmallMatrix<T,M,N,N,1,false>
  {
    public:

      inline SmallMatrixComposite() : itsm(0) {}
      inline SmallMatrixComposite(const SmallMatrixComposite<T,M,N>&) : 
	itsm(0) {}
      virtual inline ~SmallMatrixComposite() {}

      inline const T* cptr() const
      {
	if (!itsm.get()) {
	  itsm.reset(new T[M*N]);
	  AssignTom(SmallMatrixView<T,M,N,N,1,false>(itsm.get() 
	    FIRSTLAST1(itsm.get(),itsm.get()+M*N)));
	}
	return itsm.get();
      }

    private:
      mutable auto_array<T> itsm;
  };

  // These are what we want to do no matter what type Tx is:
  template <class T, size_t M, size_t N, StorageType S, IndexStyle I, class Tx> 
    inline SmallMatrix<T,M,N,S,I>& operator+=(
	SmallMatrix<T,M,N,S,I>& m, const Tx& x) 
    { m.View() += x; return m; }

  template <class T, size_t M, size_t N, StorageType S, IndexStyle I, class Tx> 
    inline SmallMatrix<T,M,N,S,I>& operator-=(
	SmallMatrix<T,M,N,S,I>& m, const Tx& x) 
    { m.View() -= x; return m; }

  template <class T, size_t M, size_t N, StorageType S, IndexStyle I, class Tx> 
    inline SmallMatrix<T,M,N,S,I>& operator*=(
	SmallMatrix<T,M,N,S,I>& m, const Tx& x) 
    { m.View() *= x; return m; }

  template <class T, size_t M, size_t N, StorageType S, IndexStyle I, class Tx> 
    inline SmallMatrix<T,M,N,S,I>& operator/=(
	SmallMatrix<T,M,N,S,I>& m, const Tx& x) 
    { m.View() /= x; return m; }

  template <class T, size_t M, size_t N, StorageType S, IndexStyle I, class Tx> 
    inline SmallMatrix<T,M,N,S,I>& operator%=(
	SmallMatrix<T,M,N,S,I>& m, const Tx& x) 
    { m.View() %= x; return m; }

  //
  // Scalar * Matrix
  //

  template <class T, class T1, size_t M, size_t N, int Si1, int Sj1, bool C1> class ProdXm : 
    public SmallMatrixComposite<T,M,N> 
  {
    public:
      inline ProdXm(const T _x, const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& _m) :
	x(_x), m(_m) {}
      inline T GetX() const { return x; }
      inline const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& GetM() const { return m; }
      inline void AssignTom(
	  const SmallMatrixView<RealType(T),M,N,N,1,false>& m0) const
      { TMVAssert(IsReal(T())); MultXM(x,m0=m); }
      inline void AssignTom(
	  const SmallMatrixView<RealType(T),M,N,1,M,false>& m0) const
      { TMVAssert(IsReal(T())); MultXM(x,m0=m); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,N,1,true>& m0) const
      { MultXM(x,m0=m); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,1,M,true>& m0) const
      { MultXM(x,m0=m); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,N,1,false>& m0) const
      { MultXM(x,m0=m); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,1,M,false>& m0) const
      { MultXM(x,m0=m); }
      inline void DoAssignTom(RealType(T)* mp,
	  const int Si, const int Sj) const
      {
	TMVAssert(IsReal(T()));
	MultXM(x,m,mp,Si,Sj);
      }
      inline void DoAssignTom(ComplexType(T)* mp,
	  const int Si, const int Sj, const bool C) const
      {
	if (C) MultXM(CONJ(x),m.Conjugate(),mp,Si,Sj);
	else MultXM(x,m,mp,Si,Sj);
      }
    private:
      const T x;
      const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& m;
  };

  // m*=x
  template <class T, size_t M, size_t N, int Si, int Sj, bool C> 
    inline const SmallMatrixView<T,M,N,Si,Sj,C>& operator*=(
	const SmallMatrixView<T,M,N,Si,Sj,C>& m, T x) 
    { MultXM(x,m); return m; }

  template <class T, size_t M, size_t N, int Si, int Sj, bool C> 
    inline const SmallMatrixView<CT,M,N,Si,Sj,C>& operator*=(
	const SmallMatrixView<CT,M,N,Si,Sj,C>& m, T x) 
    { MultXM(x,m); return m; }

  template <class T, size_t M, size_t N, int Si, int Sj, bool C> 
    inline const SmallMatrixView<CT,M,N,Si,Sj,C>& operator*=(
	const SmallMatrixView<CT,M,N,Si,Sj,C>& m, CCT x) 
    { MultXM(CT(x),m); return m; }

  template <class T, size_t M, size_t N, int Si, int Sj, bool C> 
    inline const SmallMatrixView<CT,M,N,Si,Sj,C>& operator*=(
	const SmallMatrixView<CT,M,N,Si,Sj,C>& m, VCT x) 
    { MultXM(CT(x),m); return m; }

  // m/=x
  template <class T, size_t M, size_t N, int Si, int Sj, bool C> 
    inline const SmallMatrixView<T,M,N,Si,Sj,C>& operator/=(
	const SmallMatrixView<T,M,N,Si,Sj,C>& m, T x) 
    { MultXM(T(1)/x,m); return m; }

  template <class T, size_t M, size_t N, int Si, int Sj, bool C> 
    inline const SmallMatrixView<CT,M,N,Si,Sj,C>& operator/=(
	const SmallMatrixView<CT,M,N,Si,Sj,C>& m, T x) 
    { MultXM(T(1)/x,m); return m; }

  template <class T, size_t M, size_t N, int Si, int Sj, bool C> 
    inline const SmallMatrixView<CT,M,N,Si,Sj,C>& operator/=(
	const SmallMatrixView<CT,M,N,Si,Sj,C>& m, CCT x) 
    { MultXM(T(1)/CT(x),m); return m; }

  template <class T, size_t M, size_t N, int Si, int Sj, bool C> 
    inline const SmallMatrixView<CT,M,N,Si,Sj,C>& operator/=(
	const SmallMatrixView<CT,M,N,Si,Sj,C>& m, VCT x) 
    { MultXM(T(1)/CT(x),m); return m; }

#define GENMATRIX GenSmallMatrix
#define PRODXM ProdXm
#define X ,M,N,Si1,Sj1,C1
#define Y ,size_t M, size_t N, int Si1, int Sj1, bool C1
#include "TMV_AuxProdXM.h"
  // Defines things like -m, x*m, m*x, x*(x*m), etc.
#undef PRODXM
#undef GENMATRIX

  //
  // Matrix + Scalar
  //

  template <class T, class T1, size_t N, int Si1, int Sj1, bool C1> 
    class SummX_1 : 
      public SmallMatrixComposite<T,N,N>
    {
      // m + x2
      public:
	inline SummX_1(
	    T DEBUGPARAM(_x1), const GenSmallMatrix<T1,N,N,Si1,Sj1,C1>& _m,
	    T _x2) :
	  m(_m), x2(_x2) 
	{ TMVAssert(_x1 == T(1)); }
	inline const GenSmallMatrix<T1,N,N,Si1,Sj1,C1>& GetM() const 
	{ return m; }
	inline T GetX2() const { return x2; }
	inline void AssignTom(
	    const SmallMatrixView<RealType(T),N,N,N,1,false>& m0) const
	{ TMVAssert(IsReal(T())); (m0 = m).diag().AddToAll(REAL(x2)); }
	inline void AssignTom(
	    const SmallMatrixView<RealType(T),N,N,1,N,false>& m0) const
	{ TMVAssert(IsReal(T())); (m0 = m).diag().AddToAll(REAL(x2)); }
	inline void AssignTom(
	    const SmallMatrixView<ComplexType(T),N,N,N,1,true>& m0) const
	{ (m0 = m).diag().AddToAll(ComplexType(T)(x2)); }
	inline void AssignTom(
	    const SmallMatrixView<ComplexType(T),N,N,1,N,true>& m0) const
	{ (m0 = m).diag().AddToAll(ComplexType(T)(x2)); }
	inline void AssignTom(
	    const SmallMatrixView<ComplexType(T),N,N,N,1,false>& m0) const
	{ (m0 = m).diag().AddToAll(ComplexType(T)(x2)); }
	inline void AssignTom(
	    const SmallMatrixView<ComplexType(T),N,N,1,N,false>& m0) const
	{ (m0 = m).diag().AddToAll(ComplexType(T)(x2)); }
	inline void DoAssignTom(RealType(T)* mp,
	    const int Si, const int Sj) const
	{
	  TMVAssert(IsReal(T()));
	  const int Sd = Si+Sj;
	  Copy(m,mp,Si,Sj);
	  for(size_t i=N;i>0;--i,mp+=Sd) *mp += REAL(x2);
	}
	inline void DoAssignTom(ComplexType(T)* mp,
	    const int Si, const int Sj, const bool C) const
	{
	  const int Sd = Si+Sj;
	  if (C) {
	    Copy(m.Conjugate(),mp,Si,Sj);
	    for(size_t i=N;i>0;--i,mp+=Sd) *mp += CONJ(ComplexType(T)(x2));
	  } else {
	    Copy(m,mp,Si,Sj);
	    for(size_t i=N;i>0;--i,mp+=Sd) *mp += ComplexType(T)(x2);
	  }
	}
      private:
	const GenSmallMatrix<T1,N,N,Si1,Sj1,C1>& m;
	const T x2;
    };

  template <class T, class T1, size_t N, int Si1, int Sj1, bool C1> 
    class SummX : 
      public SmallMatrixComposite<T,N,N>
    {
      // x1*m + x2
      public:
	inline SummX(T _x1, const GenSmallMatrix<T1,N,N,Si1,Sj1,C1>& _m, T _x2) :
	  x1(_x1), m(_m), x2(_x2) {}
	inline T GetX1() const { return x1; }
	inline const GenSmallMatrix<T1,N,N,Si1,Sj1,C1>& GetM() const 
	{ return m; }
	inline T GetX2() const { return x2; }
	inline void AssignTom(
	    const SmallMatrixView<RealType(T),N,N,N,1,false>& m0) const
	{ TMVAssert(IsReal(T())); (m0 = x1*m).diag().AddToAll(REAL(x2)); }
	inline void AssignTom(
	    const SmallMatrixView<RealType(T),N,N,1,N,false>& m0) const
	{ TMVAssert(IsReal(T())); (m0 = x1*m).diag().AddToAll(REAL(x2)); }
	inline void AssignTom(
	    const SmallMatrixView<ComplexType(T),N,N,N,1,true>& m0) const
	{ (m0 = x1*m).diag().AddToAll(ComplexType(T)(x2)); }
	inline void AssignTom(
	    const SmallMatrixView<ComplexType(T),N,N,1,N,true>& m0) const
	{ (m0 = x1*m).diag().AddToAll(ComplexType(T)(x2)); }
	inline void AssignTom(
	    const SmallMatrixView<ComplexType(T),N,N,N,1,false>& m0) const
	{ (m0 = x1*m).diag().AddToAll(ComplexType(T)(x2)); }
	inline void AssignTom(
	    const SmallMatrixView<ComplexType(T),N,N,1,N,false>& m0) const
	{ (m0 = x1*m).diag().AddToAll(ComplexType(T)(x2)); }
	inline void DoAssignTom(RealType(T)* mp,
	    const int Si, const int Sj) const
	{
	  TMVAssert(IsReal(T()));
	  const int Sd = Si+Sj;
	  MultXM(x1,m,mp,Si,Sj);
	  for(size_t i=N;i>0;--i,mp+=Sd) *mp += REAL(x2);
	}
	inline void DoAssignTom(ComplexType(T)* mp,
	    const int Si, const int Sj, const bool C) const
	{
	  const int Sd = Si+Sj;
	  if (C) {
	    MultXM(CONJ(x1),m.Conjugate(),mp,Si,Sj);
	    for(size_t i=N;i>0;--i,mp+=Sd) *mp += CONJ(ComplexType(T)(x2));
	  } else {
	    MultXM(x1,m,mp,Si,Sj);
	    for(size_t i=N;i>0;--i,mp+=Sd) *mp += ComplexType(T)(x2);
	  }
	}
      private:
	const T x1;
	const GenSmallMatrix<T1,N,N,Si1,Sj1,C1>& m;
	const T x2;
    };

  // m+=x
  template <class T, size_t N, int Si, int Sj, bool C> 
    inline const SmallMatrixView<T,N,N,Si,Sj,C>& operator+=(
	const SmallMatrixView<T,N,N,Si,Sj,C>& m, T x) 
    { m.diag().AddToAll(x); return m; }

  template <class T, size_t N, int Si, int Sj, bool C> 
    inline const SmallMatrixView<CT,N,N,Si,Sj,C>& operator+=(
	const SmallMatrixView<CT,N,N,Si,Sj,C>& m, T x) 
    { m.diag().AddToAll(CT(x)); return m; }

  template <class T, size_t N, int Si, int Sj, bool C> 
    inline const SmallMatrixView<CT,N,N,Si,Sj,C>& operator+=(
	const SmallMatrixView<CT,N,N,Si,Sj,C>& m, CCT x) 
    { m.diag().AddToAll(CT(x)); return m; }

  template <class T, size_t N, int Si, int Sj, bool C> 
    inline const SmallMatrixView<CT,N,N,Si,Sj,C>& operator+=(
	const SmallMatrixView<CT,N,N,Si,Sj,C>& m, VCT x) 
    { m.diag().AddToAll(CT(x)); return m; }

  // m-=x
  template <class T, size_t N, int Si, int Sj, bool C> 
    inline const SmallMatrixView<T,N,N,Si,Sj,C>& operator-=(
	const SmallMatrixView<T,N,N,Si,Sj,C>& m, T x) 
    { m.diag().AddToAll(-x); return m; }

  template <class T, size_t N, int Si, int Sj, bool C> 
    inline const SmallMatrixView<CT,N,N,Si,Sj,C>& operator-=(
	const SmallMatrixView<CT,N,N,Si,Sj,C>& m, T x) 
    { m.diag().AddToAll(CT(-x)); return m; }

  template <class T, size_t N, int Si, int Sj, bool C> 
    inline const SmallMatrixView<CT,N,N,Si,Sj,C>& operator-=(
	const SmallMatrixView<CT,N,N,Si,Sj,C>& m, CCT x) 
    { m.diag().AddToAll(-CT(x)); return m; }

  template <class T, size_t N, int Si, int Sj, bool C> 
    inline const SmallMatrixView<CT,N,N,Si,Sj,C>& operator-=(
	const SmallMatrixView<CT,N,N,Si,Sj,C>& m, VCT x) 
    { m.diag().AddToAll(-CT(x)); return m; }

#define GENMATRIX GenSmallMatrix
#define PRODXM ProdXm
#define SUMMX SummX
#define SUMMX_1 SummX_1
#define X1 ,N,N,Si1,Sj1,C1
#define X2 ,N,Si1,Sj1,C1
#define Y ,size_t N, int Si1, int Sj1, bool C1
#include "TMV_AuxSumMX.h"
  // Defines things like m+x, x+m, x-m, m-x, x+x*m, x*(x+m), etc.
#undef SUMMX
#undef GENMATRIX
#undef PRODXM

  //
  // Vector ^ Vector (OuterProduct)
  //

  template <class T, class T1, class T2, size_t M, size_t N, int S1, bool C1, int S2, bool C2> 
    class OProdvv_1 : 
      public SmallMatrixComposite<T,M,N>
    {
      public:
	inline OProdvv_1(
	    const T DEBUGPARAM(_x), const GenSmallVector<T1,M,S1,C1>& _v1,
	    const GenSmallVector<T2,N,S2,C2>& _v2) : v1(_v1), v2(_v2)
	{ TMVAssert(_x == T(1)); }
	inline const GenSmallVector<T1,M,S1,C1>& GetV1() const { return v1; }
	inline const GenSmallVector<T2,N,S2,C2>& GetV2() const { return v2; }
	inline void AssignTom(
	    const SmallMatrixView<RealType(T),M,N,N,1,false>& m0) const
	{ TMVAssert(IsReal(T())); Rank1Update_1(v1,v2,m0); }
	inline void AssignTom(
	    const SmallMatrixView<RealType(T),M,N,1,M,false>& m0) const
	{ TMVAssert(IsReal(T())); Rank1Update_1(v1,v2,m0); }
	inline void AssignTom(
	    const SmallMatrixView<ComplexType(T),M,N,N,1,true>& m0) const
	{ Rank1Update_1(v1,v2,m0); }
	inline void AssignTom(
	    const SmallMatrixView<ComplexType(T),M,N,1,M,true>& m0) const
	{ Rank1Update_1(v1,v2,m0); }
	inline void AssignTom(
	    const SmallMatrixView<ComplexType(T),M,N,N,1,false>& m0) const
	{ Rank1Update_1(v1,v2,m0); }
	inline void AssignTom(
	    const SmallMatrixView<ComplexType(T),M,N,1,M,false>& m0) const
	{ Rank1Update_1(v1,v2,m0); }
	inline void DoAssignTom(RealType(T)* mp, 
	    const int Si, const int Sj) const
	{
	  TMVAssert(IsReal(T()));
	  Rank1Update_1(v1,v2,mp,Si,Sj);
	}
	inline void DoAssignTom(ComplexType(T)* mp,
	    const int Si, const int Sj, const bool C) const
	{
	  if (C) Rank1Update_1(v1.Conjugate(),v2.Conjugate(),mp,Si,Sj);
	  else Rank1Update_1(v1,v2,mp,Si,Sj);
	}
      private:
	const GenSmallVector<T1,M,S1,C1>& v1;
	const GenSmallVector<T2,N,S2,C2>& v2;
    };

  template <class T, class T1, class T2, size_t M, size_t N, int S1, bool C1, int S2, bool C2> 
    class OProdvv : 
      public SmallMatrixComposite<T,M,N>
    {
      public:
	inline OProdvv(const T _x, const GenSmallVector<T1,M,S1,C1>& _v1,
	    const GenSmallVector<T2,N,S2,C2>& _v2) :
	  x(_x), v1(_v1), v2(_v2) {}
	inline T GetX() const { return x; }
	inline const GenSmallVector<T1,M,S1,C1>& GetV1() const { return v1; }
	inline const GenSmallVector<T2,N,S2,C2>& GetV2() const { return v2; }
	inline void AssignTom(
	    const SmallMatrixView<RealType(T),M,N,N,1,false>& m0) const
	{ TMVAssert(IsReal(T())); Rank1Update(x,v1,v2,m0); }
	inline void AssignTom(
	    const SmallMatrixView<RealType(T),M,N,1,M,false>& m0) const
	{ TMVAssert(IsReal(T())); Rank1Update(x,v1,v2,m0); }
	inline void AssignTom(
	    const SmallMatrixView<ComplexType(T),M,N,N,1,true>& m0) const
	{ Rank1Update(x,v1,v2,m0); }
	inline void AssignTom(
	    const SmallMatrixView<ComplexType(T),M,N,1,M,true>& m0) const
	{ Rank1Update(x,v1,v2,m0); }
	inline void AssignTom(
	    const SmallMatrixView<ComplexType(T),M,N,N,1,false>& m0) const
	{ Rank1Update(x,v1,v2,m0); }
	inline void AssignTom(
	    const SmallMatrixView<ComplexType(T),M,N,1,M,false>& m0) const
	{ Rank1Update(x,v1,v2,m0); }
	inline void DoAssignTom(RealType(T)* mp, 
	    const int Si, const int Sj) const
	{
	  TMVAssert(IsReal(T()));
	  Rank1Update(x,v1,v2,mp,Si,Sj);
	}
	inline void DoAssignTom(ComplexType(T)* mp,
	    const int Si, const int Sj, const bool C) const
	{
	  if (C) Rank1Update(CONJ(x),v1.Conjugate(),v2.Conjugate(),mp,Si,Sj);
	  else Rank1Update(x,v1,v2,mp,Si,Sj);
	}
      private:
	T x;
	const GenSmallVector<T1,M,S1,C1>& v1;
	const GenSmallVector<T2,N,S2,C2>& v2;
    };

  // m+=(v^v)
  template <class T, class T1, class T2, size_t M, size_t N, int Si0, int Sj0, bool C0, int S1, bool C1, int S2, bool C2> 
    inline const SmallMatrixView<T,M,N,Si0,Sj0,C0>& operator+=(
	const SmallMatrixView<T,M,N,Si0,Sj0,C0>& m0,
	const OProdvv_1<T,T1,T2,M,N,S1,C1,S2,C2>& opvv)
    { 
      TMVAssert(m0.colsize() == opvv.colsize());
      TMVAssert(m0.rowsize() == opvv.rowsize());
      AddRank1Update_1(opvv.GetV1(), opvv.GetV2(), m0);
      return m0;
    }

  template <class T, size_t M, size_t N, int Si0, int Sj0, bool C0, int S1, bool C1, int S2, bool C2> 
    inline const SmallMatrixView<CT,M,N,Si0,Sj0,C0>& operator+=(
	const SmallMatrixView<CT,M,N,Si0,Sj0,C0>& m0, 
	const OProdvv_1<T,T,T,M,N,S1,C1,S2,C2>& opvv)
    { 
      TMVAssert(m0.colsize() == opvv.colsize());
      TMVAssert(m0.rowsize() == opvv.rowsize());
      AddRank1Update_1(opvv.GetV1(), opvv.GetV2(), m0);
      return m0;
    }

  // m-=(v^v)
  template <class T, class T1, class T2, size_t M, size_t N, int Si0, int Sj0, bool C0, int S1, bool C1, int S2, bool C2> 
    inline const SmallMatrixView<T,M,N,Si0,Sj0,C0>& operator-=(
	const SmallMatrixView<T,M,N,Si0,Sj0,C0>& m0, 
	const OProdvv_1<T,T1,T2,M,N,S1,C1,S2,C2>& opvv)
    { 
      TMVAssert(m0.colsize() == opvv.colsize());
      TMVAssert(m0.rowsize() == opvv.rowsize());
      AddRank1Update_m1(opvv.GetV1(), opvv.GetV2(), m0);
      return m0;
    }

  template <class T, size_t M, size_t N, int Si0, int Sj0, bool C0, int S1, bool C1, int S2, bool C2> 
    inline const SmallMatrixView<CT,M,N,Si0,Sj0,C0>& operator-=(
	const SmallMatrixView<CT,M,N,Si0,Sj0,C0>& m0, 
	const OProdvv_1<T,T,T,M,N,S1,C1,S2,C2>& opvv)
    { 
      TMVAssert(m0.colsize() == opvv.colsize());
      TMVAssert(m0.rowsize() == opvv.rowsize());
      AddRank1Update_m1(opvv.GetV1(), opvv.GetV2(), m0);
      return m0;
    }

  // m+=(x*v^v)
  template <class T, class T1, class T2, size_t M, size_t N, int Si0, int Sj0, bool C0, int S1, bool C1, int S2, bool C2> 
    inline const SmallMatrixView<T,M,N,Si0,Sj0,C0>& operator+=(
	const SmallMatrixView<T,M,N,Si0,Sj0,C0>& m0,
	const OProdvv<T,T1,T2,M,N,S1,C1,S2,C2>& opvv)
    { 
      TMVAssert(m0.colsize() == opvv.colsize());
      TMVAssert(m0.rowsize() == opvv.rowsize());
      AddRank1Update(opvv.GetX(), opvv.GetV1(), opvv.GetV2(), m0);
      return m0;
    }

  template <class T, size_t M, size_t N, int Si0, int Sj0, bool C0, int S1, bool C1, int S2, bool C2> 
    inline const SmallMatrixView<CT,M,N,Si0,Sj0,C0>& operator+=(
	const SmallMatrixView<CT,M,N,Si0,Sj0,C0>& m0, 
	const OProdvv<T,T,T,M,N,S1,C1,S2,C2>& opvv)
    { 
      TMVAssert(m0.colsize() == opvv.colsize());
      TMVAssert(m0.rowsize() == opvv.rowsize());
      AddRank1Update(opvv.GetX(), opvv.GetV1(), opvv.GetV2(), m0);
      return m0;
    }

  template <class T, class T1, class T2, size_t M, size_t N, int Si0, int Sj0, bool C0> 
    inline const SmallMatrixView<T,M,N,Si0,Sj0,C0>& operator+=(
	const SmallMatrixView<T,M,N,Si0,Sj0,C0>& m0,
	const OProdVV<T,T1,T2>& opvv)
    { 
      TMVAssert(m0.colsize() == opvv.colsize());
      TMVAssert(m0.rowsize() == opvv.rowsize());
      Rank1Update<true>(opvv.GetX(), opvv.GetV1(), opvv.GetV2(), m0.RegView());
      return m0;
    }

  template <class T, size_t M, size_t N, int Si0, int Sj0, bool C0> 
    inline const SmallMatrixView<CT,M,N,Si0,Sj0,C0>& operator+=(
	const SmallMatrixView<CT,M,N,Si0,Sj0,C0>& m0, 
	const OProdVV<T,T,T>& opvv)
    { 
      TMVAssert(m0.colsize() == opvv.colsize());
      TMVAssert(m0.rowsize() == opvv.rowsize());
      Rank1Update<true>(opvv.GetX(), opvv.GetV1(), opvv.GetV2(), m0.RegView());
      return m0;
    }

  // m-=(x*v^v)
  template <class T, class T1, class T2, size_t M, size_t N, int Si0, int Sj0, bool C0, int S1, bool C1, int S2, bool C2> 
    inline const SmallMatrixView<T,M,N,Si0,Sj0,C0>& operator-=(
	const SmallMatrixView<T,M,N,Si0,Sj0,C0>& m0, 
	const OProdvv<T,T1,T2,M,N,S1,C1,S2,C2>& opvv)
    { 
      TMVAssert(m0.colsize() == opvv.colsize());
      TMVAssert(m0.rowsize() == opvv.rowsize());
      AddRank1Update(-opvv.GetX(), opvv.GetV1(), opvv.GetV2(), m0);
      return m0;
    }

  template <class T, size_t M, size_t N, int Si0, int Sj0, bool C0, int S1, bool C1, int S2, bool C2> 
    inline const SmallMatrixView<CT,M,N,Si0,Sj0,C0>& operator-=(
	const SmallMatrixView<CT,M,N,Si0,Sj0,C0>& m0, 
	const OProdvv<T,T,T,M,N,S1,C1,S2,C2>& opvv)
    { 
      TMVAssert(m0.colsize() == opvv.colsize());
      TMVAssert(m0.rowsize() == opvv.rowsize());
      AddRank1Update(-opvv.GetX(), opvv.GetV1(), opvv.GetV2(), m0);
      return m0;
    }

  template <class T, class T1, class T2, size_t M, size_t N, int Si0, int Sj0, bool C0> 
    inline const SmallMatrixView<T,M,N,Si0,Sj0,C0>& operator-=(
	const SmallMatrixView<T,M,N,Si0,Sj0,C0>& m0,
	const OProdVV<T,T1,T2>& opvv)
    { 
      TMVAssert(m0.colsize() == opvv.colsize());
      TMVAssert(m0.rowsize() == opvv.rowsize());
      Rank1Update<true>(-opvv.GetX(), opvv.GetV1(), opvv.GetV2(), 
	  m0.RegView());
      return m0;
    }

  template <class T, size_t M, size_t N, int Si0, int Sj0, bool C0> 
    inline const SmallMatrixView<CT,M,N,Si0,Sj0,C0>& operator-=(
	const SmallMatrixView<CT,M,N,Si0,Sj0,C0>& m0, 
	const OProdVV<T,T,T>& opvv)
    { 
      TMVAssert(m0.colsize() == opvv.colsize());
      TMVAssert(m0.rowsize() == opvv.rowsize());
      Rank1Update<true>(-opvv.GetX(), opvv.GetV1(), opvv.GetV2(), 
	  m0.RegView());
      return m0;
    }

#define OPRODVV OProdvv
#define OPRODVV_1 OProdvv_1
#define GENVECTOR1 GenSmallVector
#define GENVECTOR2 GenSmallVector
#define PRODXV1 ProdXv
#define PRODXV2 ProdXv
#define X1 ,M,S1,C1
#define X2 ,N,S2,C2
#define X3 ,M,N,S1,C1,S2,C2
#define Y ,size_t M,size_t N,int S1,bool C1,int S2,bool C2
#include "TMV_AuxOProdVV.h"
#undef OPRODVV
#undef GENVECTOR1
#undef GENVECTOR2
#undef PRODXv1
#undef PRODXv2


  //
  // Matrix + Matrix
  //

  template <class T, class T1, class T2, size_t M, size_t N, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2> class Summm_1_1 : 
    public SmallMatrixComposite<T,M,N> 
  {
    public:
      inline Summm_1_1(
	  const T DEBUGPARAM(_x1), 
	  const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& _m1, 
	  const T DEBUGPARAM(_x2),
	  const GenSmallMatrix<T2,M,N,Si2,Sj2,C2>& _m2) :
	m1(_m1), m2(_m2) 
      { TMVAssert(_x1 == T(1) && _x2 == T(1)); }
      inline const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& GetM1() const 
      { return m1; }
      inline const GenSmallMatrix<T2,M,N,Si1,Sj1,C1>& GetM2() const 
      { return m2; }
      inline void AssignTom(
	  const SmallMatrixView<RealType(T),M,N,N,1,false>& m0) const
      { TMVAssert(IsReal(T())); AddMM_1_1(m1,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<RealType(T),M,N,1,M,false>& m0) const
      { TMVAssert(IsReal(T())); AddMM_1_1(m1,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,N,1,true>& m0) const
      { AddMM_1_1(m1,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,1,M,true>& m0) const
      { AddMM_1_1(m1,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,N,1,false>& m0) const
      { AddMM_1_1(m1,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,1,M,false>& m0) const
      { AddMM_1_1(m1,m2,m0); }
      inline void DoAssignTom(RealType(T)* mp, 
	  const int Si, const int Sj) const
      { 
	TMVAssert(IsReal(T()));
	AddMM_1_1(m1,m2,mp,Si,Sj);
      }
      inline void DoAssignTom(ComplexType(T)* mp, 
	  const int Si, const int Sj, const bool C) const
      {
	if (C) AddMM_1_1(m1.Conjugate(),m2.Conjugate(),mp,Si,Sj);
	else AddMM_1_1(m1,m2,mp,Si,Sj);
      }
    private:
      const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& m1;
      const GenSmallMatrix<T2,M,N,Si2,Sj2,C2>& m2;
  };

  template <class T, class T1, class T2, size_t M, size_t N, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2> class Summm_1_m1 : 
    public SmallMatrixComposite<T,M,N> 
  {
    public:
      inline Summm_1_m1(
	  const T DEBUGPARAM(_x1), 
	  const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& _m1, 
	  const T DEBUGPARAM(_x2), 
	  const GenSmallMatrix<T2,M,N,Si2,Sj2,C2>& _m2) :
	m1(_m1), m2(_m2)
      { TMVAssert(_x1 == T(1) && _x2 == T(-1)); }
      inline const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& GetM1() const 
      { return m1; }
      inline const GenSmallMatrix<T2,M,N,Si1,Sj1,C1>& GetM2() const 
      { return m2; }
      inline void AssignTom(
	  const SmallMatrixView<RealType(T),M,N,N,1,false>& m0) const
      { TMVAssert(IsReal(T())); AddMM_1_m1(m1,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<RealType(T),M,N,1,M,false>& m0) const
      { TMVAssert(IsReal(T())); AddMM_1_m1(m1,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,N,1,true>& m0) const
      { AddMM_1_m1(m1,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,1,M,true>& m0) const
      { AddMM_1_m1(m1,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,N,1,false>& m0) const
      { AddMM_1_m1(m1,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,1,M,false>& m0) const
      { AddMM_1_m1(m1,m2,m0); }
      inline void DoAssignTom(RealType(T)* mp, 
	  const int Si, const int Sj) const
      { 
	TMVAssert(IsReal(T()));
	AddMM_1_m1(m1,m2,mp,Si,Sj);
      }
      inline void DoAssignTom(ComplexType(T)* mp, 
	  const int Si, const int Sj, const bool C) const
      {
	if (C) AddMM_1_m1(m1.Conjugate(),m2.Conjugate(),mp,Si,Sj);
	else AddMM_1_m1(m1,m2,mp,Si,Sj);
      }
    private:
      const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& m1;
      const GenSmallMatrix<T2,M,N,Si2,Sj2,C2>& m2;
  };

  template <class T, class T1, class T2, size_t M, size_t N, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2> class Summm_1_x : 
    public SmallMatrixComposite<T,M,N> 
  {
    public:
      inline Summm_1_x(
	  const T DEBUGPARAM(_x1), 
	  const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& _m1, 
	  const T _x2, const GenSmallMatrix<T2,M,N,Si2,Sj2,C2>& _m2) :
	m1(_m1), x2(_x2), m2(_m2) 
      { TMVAssert(_x1 == T(1)); }
      inline const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& GetM1() const 
      { return m1; }
      inline T GetX2() const { return x2; }
      inline const GenSmallMatrix<T2,M,N,Si1,Sj1,C1>& GetM2() const 
      { return m2; }
      inline void AssignTom(
	  const SmallMatrixView<RealType(T),M,N,N,1,false>& m0) const
      { TMVAssert(IsReal(T())); AddMM_1_x(m1,x2,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<RealType(T),M,N,1,M,false>& m0) const
      { TMVAssert(IsReal(T())); AddMM_1_x(m1,x2,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,N,1,true>& m0) const
      { AddMM_1_x(m1,x2,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,1,M,true>& m0) const
      { AddMM_1_x(m1,x2,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,N,1,false>& m0) const
      { AddMM_1_x(m1,x2,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,1,M,false>& m0) const
      { AddMM_1_x(m1,x2,m2,m0); }
      inline void DoAssignTom(RealType(T)* mp, 
	  const int Si, const int Sj) const
      { 
	TMVAssert(IsReal(T()));
	AddMM_1_x(m1,x2,m2,mp,Si,Sj);
      }
      inline void DoAssignTom(ComplexType(T)* mp, 
	  const int Si, const int Sj, const bool C) const
      {
	if (C) AddMM_1_x(m1.Conjugate(),CONJ(x2),m2.Conjugate(),mp,Si,Sj);
	else AddMM_1_x(m1,x2,m2,mp,Si,Sj);
      }
    private:
      const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& m1;
      const T x2;
      const GenSmallMatrix<T2,M,N,Si2,Sj2,C2>& m2;
  };

  template <class T, class T1, class T2, size_t M, size_t N, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2> class Summm_x_1 : 
    public SmallMatrixComposite<T,M,N> 
  {
    public:
      inline Summm_x_1(
	  const T _x1, const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& _m1, 
	  const T DEBUGPARAM(_x2), 
	  const GenSmallMatrix<T2,M,N,Si2,Sj2,C2>& _m2) :
	x1(_x1), m1(_m1), m2(_m2) 
      { TMVAssert(_x2 == T(1)); }
      inline T GetX1() const { return x1; }
      inline const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& GetM1() const 
      { return m1; }
      inline const GenSmallMatrix<T2,M,N,Si1,Sj1,C1>& GetM2() const 
      { return m2; }
      inline void AssignTom(
	  const SmallMatrixView<RealType(T),M,N,N,1,false>& m0) const
      { TMVAssert(IsReal(T())); AddMM_x_1(x1,m1,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<RealType(T),M,N,1,M,false>& m0) const
      { TMVAssert(IsReal(T())); AddMM_x_1(x1,m1,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,N,1,true>& m0) const
      { AddMM_x_1(x1,m1,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,1,M,true>& m0) const
      { AddMM_x_1(x1,m1,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,N,1,false>& m0) const
      { AddMM_x_1(x1,m1,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,1,M,false>& m0) const
      { AddMM_x_1(x1,m1,m2,m0); }
      inline void DoAssignTom(RealType(T)* mp, 
	  const int Si, const int Sj) const
      { 
	TMVAssert(IsReal(T()));
	AddMM_x_1(x1,m1,m2,mp,Si,Sj);
      }
      inline void DoAssignTom(ComplexType(T)* mp, 
	  const int Si, const int Sj, const bool C) const
      {
	if (C) AddMM_x_1(CONJ(x1),m1.Conjugate(),m2.Conjugate(),mp,Si,Sj);
	else AddMM_x_1(x1,m1,m2,mp,Si,Sj);
      }
    private:
      const T x1;
      const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& m1;
      const GenSmallMatrix<T2,M,N,Si2,Sj2,C2>& m2;
  };

  template <class T, class T1, class T2, size_t M, size_t N, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2> class Summm_x_m1 : 
    public SmallMatrixComposite<T,M,N> 
  {
    public:
      inline Summm_x_m1(
	  const T _x1, const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& _m1, 
	  const T DEBUGPARAM(_x2), 
	  const GenSmallMatrix<T2,M,N,Si2,Sj2,C2>& _m2) :
	x1(_x1), m1(_m1), m2(_m2) 
      { TMVAssert(_x2 == T(-1)); }
      inline T GetX1() const { return x1; }
      inline const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& GetM1() const 
      { return m1; }
      inline const GenSmallMatrix<T2,M,N,Si1,Sj1,C1>& GetM2() const 
      { return m2; }
      inline void AssignTom(
	  const SmallMatrixView<RealType(T),M,N,N,1,false>& m0) const
      { TMVAssert(IsReal(T())); AddMM_x_m1(x1,m1,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<RealType(T),M,N,1,M,false>& m0) const
      { TMVAssert(IsReal(T())); AddMM_x_m1(x1,m1,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,N,1,true>& m0) const
      { AddMM_x_m1(x1,m1,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,1,M,true>& m0) const
      { AddMM_x_m1(x1,m1,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,N,1,false>& m0) const
      { AddMM_x_m1(x1,m1,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,1,M,false>& m0) const
      { AddMM_x_m1(x1,m1,m2,m0); }
      inline void DoAssignTom(RealType(T)* mp, 
	  const int Si, const int Sj) const
      { 
	TMVAssert(IsReal(T()));
	AddMM_x_m1(x1,m1,m2,mp,Si,Sj);
      }
      inline void DoAssignTom(ComplexType(T)* mp, 
	  const int Si, const int Sj, const bool C) const
      {
	if (C) AddMM_x_m1(CONJ(x1),m1.Conjugate(),m2.Conjugate(),mp,Si,Sj);
	else AddMM_x_m1(x1,m1,m2,mp,Si,Sj);
      }
    private:
      const T x1;
      const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& m1;
      const GenSmallMatrix<T2,M,N,Si2,Sj2,C2>& m2;
  };

  template <class T, class T1, class T2, size_t M, size_t N, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2> class Summm : 
    public SmallMatrixComposite<T,M,N> 
  {
    public:
      inline Summm(const T _x1, const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& _m1, 
	  const T _x2, const GenSmallMatrix<T2,M,N,Si2,Sj2,C2>& _m2) :
	x1(_x1), m1(_m1), x2(_x2), m2(_m2) {}
      inline T GetX1() const { return x1; }
      inline const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& GetM1() const 
      { return m1; }
      inline T GetX2() const { return x2; }
      inline const GenSmallMatrix<T2,M,N,Si1,Sj1,C1>& GetM2() const 
      { return m2; }
      inline void AssignTom(
	  const SmallMatrixView<RealType(T),M,N,N,1,false>& m0) const
      { TMVAssert(IsReal(T())); AddMM(x1,m1,x2,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<RealType(T),M,N,1,M,false>& m0) const
      { TMVAssert(IsReal(T())); AddMM(x1,m1,x2,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,N,1,true>& m0) const
      { AddMM(x1,m1,x2,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,1,M,true>& m0) const
      { AddMM(x1,m1,x2,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,N,1,false>& m0) const
      { AddMM(x1,m1,x2,m2,m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,1,M,false>& m0) const
      { AddMM(x1,m1,x2,m2,m0); }
      inline void DoAssignTom(RealType(T)* mp, 
	  const int Si, const int Sj) const
      { 
	TMVAssert(IsReal(T()));
	AddMM(x1,m1,x2,m2,mp,Si,Sj);
      }
      inline void DoAssignTom(ComplexType(T)* mp, 
	  const int Si, const int Sj, const bool C) const
      {
	if (C) AddMM(CONJ(x1),m1.Conjugate(),CONJ(x2),m2.Conjugate(),mp,Si,Sj);
	else AddMM(x1,m1,x2,m2,mp,Si,Sj);
      }
    private:
      const T x1;
      const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& m1;
      const T x2;
      const GenSmallMatrix<T2,M,N,Si2,Sj2,C2>& m2;
  };

  template <class T, size_t M, size_t N, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2> 
    inline const SmallMatrixView<T,M,N,Si1,Sj1,C1>& operator+=(
	const SmallMatrixView<T,M,N,Si1,Sj1,C1>& m1, 
	const GenSmallMatrix<T,M,N,Si2,Sj2,C2>& m2) 
    { AddMM_1(m2,m1); return m1; }

  template <class T, size_t M, size_t N, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2> 
    inline const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& operator+=(
	const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& m1, 
	const GenSmallMatrix<T,M,N,Si2,Sj2,C2>& m2) 
    { AddMM_1(m2,m1); return m1; }

  template <class T, size_t M, size_t N, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2> 
    inline const SmallMatrixView<T,M,N,Si1,Sj1,C1>& operator-=(
	const SmallMatrixView<T,M,N,Si1,Sj1,C1>& m1,
	const GenSmallMatrix<T,M,N,Si2,Sj2,C2>& m2) 
    { AddMM_m1(m2,m1); return m1; }

  template <class T, size_t M, size_t N, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2> 
    inline const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& operator-=(
	const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& m1,
	const GenSmallMatrix<T,M,N,Si2,Sj2,C2>& m2) 
    { AddMM_m1(m2,m1); return m1; }

  template <class T, class T2, size_t M, size_t N, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2> 
    inline const SmallMatrixView<T,M,N,Si1,Sj1,C1>& operator+=(
	const SmallMatrixView<T,M,N,Si1,Sj1,C1>& m,
	const ProdXm<T,T2,M,N,Si2,Sj2,C2>& pxm)
    { AddMM(pxm.GetX(),pxm.GetM(),m); return m; }

  template <class T, size_t M, size_t N, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2> 
    inline const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& operator+=(
	const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& m,
	const ProdXm<T,T,M,N,Si2,Sj2,C2>& pxm)
    { AddMM(pxm.GetX(),pxm.GetM(),m); return m; }

  template <class T, class T2, size_t M, size_t N, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2> 
    inline const SmallMatrixView<T,M,N,Si1,Sj1,C1>& operator-=(
	const SmallMatrixView<T,M,N,Si1,Sj1,C1>& m,
	const ProdXm<T,T2,M,N,Si2,Sj2,C2>& pxm)
    { AddMM(-pxm.GetX(),pxm.GetM(),m); return m; }

  template <class T, size_t M, size_t N, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2> 
    inline const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& operator-=(
	const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& m,
	const ProdXm<T,T,M,N,Si2,Sj2,C2>& pxm)
    { AddMM(-pxm.GetX(),pxm.GetM(),m); return m; }

  
  template <class T, size_t M, size_t N, int Si1, int Sj1, bool C1> 
    inline const SmallMatrixView<T,M,N,Si1,Sj1,C1>& operator+=(
	const SmallMatrixView<T,M,N,Si1,Sj1,C1>& m1, const GenMatrix<T>& m2) 
    { 
      TMVAssert(m2.colsize()==M);
      TMVAssert(m2.rowsize()==N);
      AddMM(T(1),m2,m1.RegView()); 
      return m1; 
    }

  template <class T, size_t M, size_t N, int Si1, int Sj1, bool C1> 
    inline const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& operator+=(
	const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& m1, 
	const GenMatrix<T>& m2) 
    {
      TMVAssert(m2.colsize()==M);
      TMVAssert(m2.rowsize()==N);
      AddMM(T(1),m2,m1.RegView()); 
      return m1; 
    }

  template <class T, size_t M, size_t N, int Si1, int Sj1, bool C1> 
    inline const SmallMatrixView<T,M,N,Si1,Sj1,C1>& operator-=(
	const SmallMatrixView<T,M,N,Si1,Sj1,C1>& m1, const GenMatrix<T>& m2) 
    { 
      TMVAssert(m2.colsize()==M);
      TMVAssert(m2.rowsize()==N);
      AddMM(T(-1),m2,m1.RegView()); 
      return m1; 
    }

  template <class T, size_t M, size_t N, int Si1, int Sj1, bool C1> 
    inline const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& operator-=(
	const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& m1, const GenMatrix<T>& m2) 
    { 
      TMVAssert(m2.colsize()==M);
      TMVAssert(m2.rowsize()==N);
      AddMM(T(-1),m2,m1.RegView()); 
      return m1; 
    }

  template <class T, class T2, size_t M, size_t N, int Si1, int Sj1, bool C1> 
    inline const SmallMatrixView<T,M,N,Si1,Sj1,C1>& operator+=(
	const SmallMatrixView<T,M,N,Si1,Sj1,C1>& m, const ProdXM<T,T2>& pxm)
    { 
      TMVAssert(pxm.colsize()==M);
      TMVAssert(pxm.rowsize()==N);
      AddMM(pxm.GetX(),pxm.GetM(),m.RegView()); 
      return m; 
    }

  template <class T, size_t M, size_t N, int Si1, int Sj1, bool C1> 
    inline const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& operator+=(
	const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& m, const ProdXM<T,T>& pxm)
    { 
      TMVAssert(pxm.colsize()==M);
      TMVAssert(pxm.rowsize()==N);
      AddMM(pxm.GetX(),pxm.GetM(),m.RegView()); 
      return m; 
    }

  template <class T, class T2, size_t M, size_t N, int Si1, int Sj1, bool C1> 
    inline const SmallMatrixView<T,M,N,Si1,Sj1,C1>& operator-=(
	const SmallMatrixView<T,M,N,Si1,Sj1,C1>& m, const ProdXM<T,T2>& pxm)
    { 
      TMVAssert(pxm.colsize()==M);
      TMVAssert(pxm.rowsize()==N);
      AddMM(-pxm.GetX(),pxm.GetM(),m.RegView()); 
      return m; 
    }

  template <class T, size_t M, size_t N, int Si1, int Sj1, bool C1> 
    inline const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& operator-=(
	const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& m, const ProdXM<T,T>& pxm)
    { 
      TMVAssert(pxm.colsize()==M);
      TMVAssert(pxm.rowsize()==N);
      AddMM(-pxm.GetX(),pxm.GetM(),m.RegView()); 
      return m; 
    }

  
#define GENMATRIX1 GenSmallMatrix
#define GENMATRIX2 GenSmallMatrix
#define PRODXM1 ProdXm
#define PRODXM2 ProdXm
#define SUMMM Summm
#define SUMMM_1_1 Summm_1_1
#define SUMMM_1_m1 Summm_1_m1
#define SUMMM_1_x Summm_1_x
#define SUMMM_x_1 Summm_x_1
#define SUMMM_x_m1 Summm_x_m1
#define X1 ,M,N,Si1,Sj1,C1
#define X2 ,M,N,Si2,Sj2,C2
#define X3 ,M,N,Si1,Sj1,C1,Si2,Sj2,C2
#define Y ,size_t M,size_t N,int Si1,int Sj1,bool C1,int Si2,int Sj2,bool C2
#include "TMV_AuxSumMM.h"
  // Defines things like m+m, m-m, (x*m)-(x*m), etc.
  
#define SUMMM_1_1 Summm_1_1
#define SUMMM_1_m1 Summm_1_m1
#define SUMMM_1_x Summm_1_x
#define SUMMM_x_1 Summm_x_1
#define SUMMM_x_m1 Summm_x_m1
#define X3 ,M,N,Si1,Sj1,C1,Si2,Sj2,C2
#define Y ,size_t M,size_t N,int Si1,int Sj1,bool C1,int Si2,int Sj2,bool C2
#include "TMV_AuxSumMMa.h"
  // Defines things like -(m+m), x*(m+m), (m+m)/x, etc.
#undef GENMATRIX1
#undef GENMATRIX2
#undef SUMMM
#undef PRODXM1
#undef PRODXM2

  //
  // Matrix * Matrix
  //

  template <class T, class T1, class T2, size_t M, size_t N, size_t K, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2> class Prodmm_1 : 
    public SmallMatrixComposite<T,M,N>
  {
    public:
      inline Prodmm_1(
	  const T DEBUGPARAM(_x), const GenSmallMatrix<T1,M,K,Si1,Sj1,C1>& _m1, 
	  const GenSmallMatrix<T2,K,N,Si2,Sj2,C2>& _m2) :
	m1(_m1), m2(_m2) 
      { TMVAssert(_x == T(1)); }
      inline const GenSmallMatrix<T1,M,K,Si1,Sj1,C1>& GetM1() const 
      { return m1; }
      inline const GenSmallMatrix<T2,K,N,Si2,Sj2,C2>& GetM2() const 
      { return m2; }
      inline void AssignTom(
	  const SmallMatrixView<RealType(T),M,N,1,M,false>& m0) const
      { TMVAssert(IsReal(T())); MultMM_1(m1, m2, m0); }
      inline void AssignTom(
	  const SmallMatrixView<RealType(T),M,N,N,1,false>& m0) const
      { TMVAssert(IsReal(T())); MultMM_1(m1, m2, m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,1,M,true>& m0) const
      { MultMM_1(m1, m2, m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,N,1,true>& m0) const
      { MultMM_1(m1, m2, m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,1,M,false>& m0) const
      { MultMM_1(m1, m2, m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,N,1,false>& m0) const
      { MultMM_1(m1, m2, m0); }
      inline void DoAssignTom(RealType(T)* mp,
	  const int Si, const int Sj) const
      {
	TMVAssert(IsReal(T()));
	MultMM_1(m1,m2,mp,Si,Sj); 
      }
      inline void DoAssignTom(ComplexType(T)* mp,
	  const int Si, const int Sj, const bool C) const
      {
	if (C) MultMM_1(m1.Conjugate(),m2.Conjugate(),mp,Si,Sj);
	else MultMM_1(m1,m2,mp,Si,Sj); 
      }
    private:
      const GenSmallMatrix<T1,M,K,Si1,Sj1,C1>& m1;
      const GenSmallMatrix<T2,K,N,Si2,Sj2,C2>& m2;
  };

  template <class T, class T1, class T2, size_t M, size_t N, size_t K, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2> class Prodmm : 
    public SmallMatrixComposite<T,M,N>
  {
    public:
      inline Prodmm(const T _x, const GenSmallMatrix<T1,M,K,Si1,Sj1,C1>& _m1, 
	  const GenSmallMatrix<T2,K,N,Si2,Sj2,C2>& _m2) :
	x(_x), m1(_m1), m2(_m2) {}
      inline T GetX() const { return x; }
      inline const GenSmallMatrix<T1,M,K,Si1,Sj1,C1>& GetM1() const 
      { return m1; }
      inline const GenSmallMatrix<T2,K,N,Si2,Sj2,C2>& GetM2() const 
      { return m2; }
      inline void AssignTom(
	  const SmallMatrixView<RealType(T),M,N,1,M,false>& m0) const
      { TMVAssert(IsReal(T())); MultMM(x, m1, m2, m0); }
      inline void AssignTom(
	  const SmallMatrixView<RealType(T),M,N,N,1,false>& m0) const
      { TMVAssert(IsReal(T())); MultMM(x, m1, m2, m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,1,M,true>& m0) const
      { MultMM(x, m1, m2, m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,N,1,true>& m0) const
      { MultMM(x, m1, m2, m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,1,M,false>& m0) const
      { MultMM(x, m1, m2, m0); }
      inline void AssignTom(
	  const SmallMatrixView<ComplexType(T),M,N,N,1,false>& m0) const
      { MultMM(x, m1, m2, m0); }
      inline void DoAssignTom(RealType(T)* mp,
	  const int Si, const int Sj) const
      {
	TMVAssert(IsReal(T()));
	MultMM(x,m1,m2,mp,Si,Sj); 
      }
      inline void DoAssignTom(ComplexType(T)* mp,
	  const int Si, const int Sj, const bool C) const
      {
	if (C) MultMM(CONJ(x),m1.Conjugate(),m2.Conjugate(),mp,Si,Sj);
	else MultMM(x,m1,m2,mp,Si,Sj); 
      }
    private:
      const T x;
      const GenSmallMatrix<T1,M,K,Si1,Sj1,C1>& m1;
      const GenSmallMatrix<T2,K,N,Si2,Sj2,C2>& m2;
  };

  template <class T, size_t M, size_t N, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2>
    inline const SmallMatrixView<T,M,N,Si1,Sj1,C1>& operator*=(
	const SmallMatrixView<T,M,N,Si1,Sj1,C1>& m1, 
	const GenSmallMatrix<T,M,N,Si2,Sj2,C2>& m2)
    { MultMM_1(SmallMatrix<T,M,N,RowMajor>(m1),m2,m1); return m1; }

  template <class T, size_t M, size_t N, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2>
    inline const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& operator*=(
	const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& m1, 
	const GenSmallMatrix<T,M,N,Si2,Sj2,C2>& m2)
    { MultMM_1(SmallMatrix<CT,M,N,RowMajor>(m1),m2,m1); return m1; }

  template <class T, class T2, size_t M, size_t N, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2>
    inline const SmallMatrixView<T,M,N,Si1,Sj1,C1>& operator*=(
	const SmallMatrixView<T,M,N,Si1,Sj1,C1>& m1, 
	const ProdXm<T,T2,M,N,Si2,Sj2,C2>& pxm)
    { 
      MultMM(pxm.GetX(),SmallMatrix<T,M,N,RowMajor>(m1),pxm.GetM(),m1); 
      return m1; 
    }

  template <class T, size_t M, size_t N, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2>
    inline const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& operator*=(
	const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& m1, 
	const ProdXm<T,T,M,N,Si2,Sj2,C2>& pxm)
    { 
      MultMM(pxm.GetX(),SmallMatrix<CT,M,N,RowMajor>(m1),pxm.GetM(),m1); 
      return m1; 
    }

  template <class T, class T2, class T3, size_t M, size_t N, size_t K, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2, int Si3, int Sj3, bool C3>
    inline const SmallMatrixView<T,M,N,Si1,Sj1,C1>& operator+=(
	const SmallMatrixView<T,M,N,Si1,Sj1,C1>& m, 
	const Prodmm_1<T,T2,T3,M,N,K,Si2,Sj2,C2,Si3,Sj3,C3>& pmm)
    { AddMultMM_1(pmm.GetM1(),pmm.GetM2(),m); return m; }

  template <class T, size_t M, size_t N, size_t K, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2, int Si3, int Sj3, bool C3>
    inline const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& operator+=(
	const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& m, 
	const Prodmm_1<T,T,T,M,N,K,Si2,Sj2,C2,Si3,Sj3,C3>& pmm)
    { AddMultMM_1(pmm.GetM1(),pmm.GetM2(),m); return m; }

  template <class T, class T2, class T3, size_t M, size_t N, size_t K, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2, int Si3, int Sj3, bool C3>
    inline const SmallMatrixView<T,M,N,Si1,Sj1,C1>& operator-=(
	const SmallMatrixView<T,M,N,Si1,Sj1,C1>& m, 
	const Prodmm_1<T,T2,T3,M,N,K,Si2,Sj2,C2,Si3,Sj3,C3>& pmm)
    { AddMultMM_m1(pmm.GetM1(),pmm.GetM2(),m); return m; }

  template <class T, size_t M, size_t N, size_t K, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2, int Si3, int Sj3, bool C3>
    inline const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& operator-=(
	const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& m, 
	const Prodmm_1<T,T,T,M,N,K,Si2,Sj2,C2,Si3,Sj3,C3>& pmm)
    { AddMultMM_m1(pmm.GetM1(),pmm.GetM2(),m); return m; }

  template <class T, class T2, class T3, size_t M, size_t N, size_t K, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2, int Si3, int Sj3, bool C3>
    inline const SmallMatrixView<T,M,N,Si1,Sj1,C1>& operator+=(
	const SmallMatrixView<T,M,N,Si1,Sj1,C1>& m, 
	const Prodmm<T,T2,T3,M,N,K,Si2,Sj2,C2,Si3,Sj3,C3>& pmm)
    { AddMultMM(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); return m; }

  template <class T, size_t M, size_t N, size_t K, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2, int Si3, int Sj3, bool C3>
    inline const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& operator+=(
	const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& m, 
	const Prodmm<T,T,T,M,N,K,Si2,Sj2,C2,Si3,Sj3,C3>& pmm)
    { AddMultMM(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); return m; }

  template <class T, class T2, class T3, size_t M, size_t N, size_t K, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2, int Si3, int Sj3, bool C3>
    inline const SmallMatrixView<T,M,N,Si1,Sj1,C1>& operator-=(
	const SmallMatrixView<T,M,N,Si1,Sj1,C1>& m, 
	const Prodmm<T,T2,T3,M,N,K,Si2,Sj2,C2,Si3,Sj3,C3>& pmm)
    { AddMultMM(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); return m; }

  template <class T, size_t M, size_t N, size_t K, int Si1, int Sj1, bool C1, int Si2, int Sj2, bool C2, int Si3, int Sj3, bool C3>
    inline const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& operator-=(
	const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& m, 
	const Prodmm<T,T,T,M,N,K,Si2,Sj2,C2,Si3,Sj3,C3>& pmm)
    { AddMultMM(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); return m; }

  template <class T, size_t M, size_t N, int Si1, int Sj1, bool C1>
    inline const SmallMatrixView<T,M,N,Si1,Sj1,C1>& operator*=(
	const SmallMatrixView<T,M,N,Si1,Sj1,C1>& m1, 
	const GenMatrix<T>& m2)
    { 
      TMVAssert(m2.colsize() == N);
      TMVAssert(m2.rowsize() == N);
      MultMM<false>(T(1),m1,m2,m1.RegView()); 
      return m1; 
    }

  template <class T, size_t M, size_t N, int Si1, int Sj1, bool C1>
    inline const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& operator*=(
	const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& m1, 
	const GenMatrix<T>& m2)
    {
      TMVAssert(m2.colsize() == N);
      TMVAssert(m2.rowsize() == N);
      MultMM<false>(T(1),m1,m2,m1.RegView()); 
      return m1; 
    }

  template <class T, class T2, size_t M, size_t N, int Si1, int Sj1, bool C1>
    inline const SmallMatrixView<T,M,N,Si1,Sj1,C1>& operator*=(
	const SmallMatrixView<T,M,N,Si1,Sj1,C1>& m1, 
	const ProdXM<T,T2>& pxm)
    { 
      TMVAssert(pxm.colsize() == N);
      TMVAssert(pxm.rowsize() == N);
      MultMM<false>(pxm.GetX(),m1,pxm.GetM(),m1.RegView()); 
      return m1; 
    }

  template <class T, size_t M, size_t N, int Si1, int Sj1, bool C1>
    inline const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& operator*=(
	const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& m1, 
	const ProdXM<T,T>& pxm)
    { 
      TMVAssert(pxm.colsize() == N);
      TMVAssert(pxm.rowsize() == N);
      MultMM<false>(pxm.GetX(),m1,pxm.GetM(),m1.RegView()); 
      return m1; 
    }

  template <class T, class T2, class T3, size_t M, size_t N, size_t K, int Si1, int Sj1, bool C1>
    inline const SmallMatrixView<T,M,N,Si1,Sj1,C1>& operator+=(
	const SmallMatrixView<T,M,N,Si1,Sj1,C1>& m, 
	const ProdMM<T,T2,T3>& pmm)
    {
      TMVAssert(pmm.colsize() == M);
      TMVAssert(pmm.rowsize() == N);
      MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m.RegView()); 
      return m; 
    }

  template <class T, size_t M, size_t N, size_t K, int Si1, int Sj1, bool C1>
    inline const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& operator+=(
	const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& m, 
	const ProdMM<T,T,T>& pmm)
    {
      TMVAssert(pmm.colsize() == M);
      TMVAssert(pmm.rowsize() == N);
      MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m.RegView()); 
      return m; 
    }

  template <class T, class T2, class T3, size_t M, size_t N, size_t K, int Si1, int Sj1, bool C1>
    inline const SmallMatrixView<T,M,N,Si1,Sj1,C1>& operator-=(
	const SmallMatrixView<T,M,N,Si1,Sj1,C1>& m, 
	const ProdMM<T,T2,T3>& pmm)
    { 
      TMVAssert(pmm.colsize() == M);
      TMVAssert(pmm.rowsize() == N);
      MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m.RegView()); 
      return m; 
    }

  template <class T, size_t M, size_t N, size_t K, int Si1, int Sj1, bool C1>
    inline const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& operator-=(
	const SmallMatrixView<CT,M,N,Si1,Sj1,C1>& m, 
	const ProdMM<T,T,T>& pmm)
    {
      TMVAssert(pmm.colsize() == M);
      TMVAssert(pmm.rowsize() == N);
      MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m.RegView()); 
      return m; 
    }

#define GENMATRIX1 GenSmallMatrix
#define GENMATRIX2 GenSmallMatrix
#define PRODMM Prodmm
#define PRODMM_1 Prodmm_1
#define PRODXM1 ProdXm
#define PRODXM2 ProdXm
#define X1 ,M,K,Si1,Sj1,C1
#define X2 ,K,N,Si2,Sj2,C2
#define X3 ,M,N,K,Si1,Sj1,C1,Si2,Sj2,C2
#define Y ,size_t M,size_t N,size_t K,int Si1,int Sj1,bool C1,int Si2,int Sj2,bool C2
#include "TMV_AuxProdMM.h"
  // Defines things like m*m, m*(x*m), etc.
#define PRODMM_1 Prodmm_1
#define X3 ,M,N,K,Si1,Sj1,C1,Si2,Sj2,C2
#define Y ,size_t M,size_t N,size_t K,int Si1,int Sj1,bool C1,int Si2,int Sj2,bool C2
#include "TMV_AuxProdMMa.h"
  // Defines things like -(m*m), x*(m*m), etc.
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2


  //
  // Matrix * Vector
  //

  template <class T, class T1, class T2, size_t M, size_t N, int Si1, int Sj1, bool C1, int S2, bool C2> class Prodmv_1 :
    public SmallVectorComposite<T,M>
  {
    public:
      inline Prodmv_1(
	  const T DEBUGPARAM(_x), const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& _m,
	  const GenSmallVector<T2,N,S2,C2>& _v) :
	m(_m), v(_v) 
      { TMVAssert(_x == T(1)); }
      inline const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& GetM() const { return m; }
      inline const GenSmallVector<T2,N,S2,C2>& GetV() const { return v; }
      inline void AssignTov(
	  const SmallVectorView<RealType(T),M,1,false>& v0) const
      { TMVAssert(IsReal(T())); MultMV_1(m,v,v0); }
      inline void AssignTov(
	  const SmallVectorView<ComplexType(T),M,1,true>& v0) const
      { MultMV_1(m,v,v0); }
      inline void AssignTov(
	  const SmallVectorView<ComplexType(T),M,1,false>& v0) const
      { MultMV_1(m,v,v0); }
      inline void DoAssignTov(RealType(T)* vp, const int S) const
      { TMVAssert(IsReal(T())); MultMV_1(m,v,vp,S); }
      inline void DoAssignTov(ComplexType(T)* vp,
	  const int S, const bool C) const
      { 
	if (C) MultMV_1(m.Conjugate(),v.Conjugate(),vp,S); 
	else MultMV_1(m,v,vp,S); 
      }
    private:
      const ConstSmallMatrixView<T1,M,N,Si1,Sj1,C1> m;
      const GenSmallVector<T2,N,S2,C2>& v;
  };

  template <class T, class T1, class T2, size_t M, size_t N, int Si1, int Sj1, bool C1, int S2, bool C2> class Prodmv :
    public SmallVectorComposite<T,M>
  {
    public:
      inline Prodmv(const T _x, const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& _m,
	  const GenSmallVector<T2,N,S2,C2>& _v) :
	x(_x), m(_m), v(_v) {}
      inline T GetX() const { return x; }
      inline const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& GetM() const { return m; }
      inline const GenSmallVector<T2,N,S2,C2>& GetV() const { return v; }
      inline void AssignTov(
	  const SmallVectorView<RealType(T),M,1,false>& v0) const
      { TMVAssert(IsReal(T())); MultMV(x,m,v,v0); }
      inline void AssignTov(
	  const SmallVectorView<ComplexType(T),M,1,true>& v0) const
      { MultMV(x,m,v,v0); }
      inline void AssignTov(
	  const SmallVectorView<ComplexType(T),M,1,false>& v0) const
      { MultMV(x,m,v,v0); }
      inline void DoAssignTov(RealType(T)* vp, const int S) const
      { TMVAssert(IsReal(T())); MultMV(x,m,v,vp,S); }
      inline void DoAssignTov(ComplexType(T)* vp,
	  const int S, const bool C) const
      { 
	if (C) MultMV(CONJ(x),m.Conjugate(),v.Conjugate(),vp,S); 
	else MultMV(x,m,v,vp,S); 
      }
    private:
      const T x;
      const ConstSmallMatrixView<T1,M,N,Si1,Sj1,C1> m;
      const GenSmallVector<T2,N,S2,C2>& v;
  };

  template <class T, size_t N, int Si1, int Sj1, bool C1, int S2, bool C2> 
    inline const SmallVectorView<T,N,S2,C2>& operator*=(
	const SmallVectorView<T,N,S2,C2>& v,
	const GenSmallMatrix<T,N,N,Si1,Sj1,C1>& m)
    { MultMV_1(m.Transpose(),SmallVector<T,N>(v),v); return v; }

  template <class T, size_t N, int Si1, int Sj1, bool C1, int S2, bool C2> 
    inline const SmallVectorView<CT,N,S2,C2>& operator*=(
	const SmallVectorView<CT,N,S2,C2>& v,
	const GenSmallMatrix<T,N,N,Si1,Sj1,C1>& m)
    { MultMV_1(m.Transpose(),SmallVector<CT,N>(v),v); return v; }

  template <class T, class T1, size_t N, int Si1, int Sj1, bool C1, int S2, bool C2> 
    inline const SmallVectorView<T,N,S2,C2>& operator*=(
	const SmallVectorView<T,N,S2,C2>& v,
	const ProdXm<T,T1,N,N,Si1,Sj1,C1>& pxm)
    { 
      MultMV(pxm.GetX(),pxm.GetM().Transpose(),SmallVector<T,N>(v),v);
      return v;
    }

  template <class T, size_t N, int Si1, int Sj1, bool C1, int S2, bool C2> 
    inline const SmallVectorView<CT,N,S2,C2>& operator*=(
	const SmallVectorView<CT,N,S2,C2>& v,
	const ProdXm<T,T,N,N,Si1,Sj1,C1>& pxm)
    {
      MultMV(pxm.GetX(),pxm.GetM().Transpose(),SmallVector<CT,N>(v),v);
      return v;
    }

  template <class T, class T1, class T2, size_t M, size_t N, int Si1, int Sj1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline const SmallVectorView<T,M,S3,C3>& operator+=(
	const SmallVectorView<T,M,S3,C3>& v, 
	const Prodmv_1<T,T1,T2,M,N,Si1,Sj1,C1,S2,C2>& pmv)
    { AddMultMV_1(pmv.GetM(),pmv.GetV(),v); return v; }

  template <class T, size_t M, size_t N, int Si1, int Sj1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline const SmallVectorView<CT,M,S3,C3>& operator+=(
	const SmallVectorView<CT,M,S3,C3>& v,
	const Prodmv_1<T,T,T,M,N,Si1,Sj1,C1,S2,C2>& pmv)
    { AddMultMV_1(pmv.GetM(),pmv.GetV(),v); return v; }

  template <class T, class T1, class T2, size_t M, size_t N, int Si1, int Sj1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline const SmallVectorView<T,M,S3,C3>& operator-=(
	const SmallVectorView<T,M,S3,C3>& v,
	const Prodmv_1<T,T1,T2,M,N,Si1,Sj1,C1,S2,C2>& pmv)
    { AddMultMV_m1(pmv.GetM(), pmv.GetV(), v); return v; }

  template <class T, size_t M, size_t N, int Si1, int Sj1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline const SmallVectorView<CT,M,S3,C3>& operator-=(
	const SmallVectorView<CT,M,S3,C3>& v,
	const Prodmv_1<T,T,T,M,N,Si1,Sj1,C1,S2,C2>& pmv)
    { AddMultMV_m1(pmv.GetM(),pmv.GetV(),v); return v; }

  template <class T, class T1, class T2, size_t M, size_t N, int Si1, int Sj1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline const SmallVectorView<T,M,S3,C3>& operator+=(
	const SmallVectorView<T,M,S3,C3>& v, 
	const Prodmv<T,T1,T2,M,N,Si1,Sj1,C1,S2,C2>& pmv)
    { AddMultMV(pmv.GetX(),pmv.GetM(),pmv.GetV(),v); return v; }

  template <class T, size_t M, size_t N, int Si1, int Sj1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline const SmallVectorView<CT,M,S3,C3>& operator+=(
	const SmallVectorView<CT,M,S3,C3>& v,
	const Prodmv<T,T,T,M,N,Si1,Sj1,C1,S2,C2>& pmv)
    { AddMultMV(pmv.GetX(),pmv.GetM(),pmv.GetV(),v); return v; }

  template <class T, class T1, class T2, size_t M, size_t N, int Si1, int Sj1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline const SmallVectorView<T,M,S3,C3>& operator-=(
	const SmallVectorView<T,M,S3,C3>& v,
	const Prodmv<T,T1,T2,M,N,Si1,Sj1,C1,S2,C2>& pmv)
    { AddMultMV(-pmv.GetX(), pmv.GetM(), pmv.GetV(), v); return v; }

  template <class T, size_t M, size_t N, int Si1, int Sj1, bool C1, int S2, bool C2, int S3, bool C3> 
    inline const SmallVectorView<CT,M,S3,C3>& operator-=(
	const SmallVectorView<CT,M,S3,C3>& v,
	const Prodmv<T,T,T,M,N,Si1,Sj1,C1,S2,C2>& pmv)
    { AddMultMV(-pmv.GetX(),pmv.GetM(),pmv.GetV(),v); return v; }

  template <class T, size_t N, int S2, bool C2> 
    inline const SmallVectorView<T,N,S2,C2>& operator*=(
	const SmallVectorView<T,N,S2,C2>& v, const GenMatrix<T>& m)
    { 
      TMVAssert(m.colsize() == N);
      TMVAssert(m.rowsize() == N);
      MultMV<false>(T(1),m.Transpose(),SmallVector<T,N>(v),v); 
      return v; 
    }

  template <class T, size_t N, int S2, bool C2> 
    inline const SmallVectorView<CT,N,S2,C2>& operator*=(
	const SmallVectorView<CT,N,S2,C2>& v, const GenMatrix<T>& m)
    {
      TMVAssert(m.colsize() == N);
      TMVAssert(m.rowsize() == N);
      MultMV<false>(T(1),m.Transpose(),SmallVector<CT,N>(v),v); 
      return v; 
    }

  template <class T, class T1, size_t N, int S2, bool C2> 
    inline const SmallVectorView<T,N,S2,C2>& operator*=(
	const SmallVectorView<T,N,S2,C2>& v, const ProdXM<T,T1>& pxm)
    { 
      TMVAssert(pxm.colsize() == N);
      TMVAssert(pxm.rowsize() == N);
      MultMV<false>(pxm.GetX(),pxm.GetM().Transpose(),SmallVector<T,N>(v),v);
      return v;
    }

  template <class T, size_t N, int S2, bool C2> 
    inline const SmallVectorView<CT,N,S2,C2>& operator*=(
	const SmallVectorView<CT,N,S2,C2>& v, const ProdXM<T,T>& pxm)
    {
      TMVAssert(pxm.colsize() == N);
      TMVAssert(pxm.rowsize() == N);
      MultMV<false>(pxm.GetX(),pxm.GetM().Transpose(),SmallVector<CT,N>(v),v);
      return v;
    }

  template <class T, class T1, class T2, size_t M, size_t N, int S3, bool C3> 
    inline const SmallVectorView<T,M,S3,C3>& operator+=(
	const SmallVectorView<T,M,S3,C3>& v, 
	const ProdMV<T,T1,T2>& pmv)
    {
      TMVAssert(pmv.size() == N);
      MultMV<true>(pmv.GetX(),pmv.GetM(),pmv.GetV(),v); 
      return v; 
    }

  template <class T, size_t M, size_t N, int S3, bool C3> 
    inline const SmallVectorView<CT,M,S3,C3>& operator+=(
	const SmallVectorView<CT,M,S3,C3>& v,
	const ProdMV<T,T,T>& pmv)
    {
      TMVAssert(pmv.size() == N);
      MultMV<true>(pmv.GetX(),pmv.GetM(),pmv.GetV(),v); 
      return v; 
    }

  template <class T, class T1, class T2, size_t M, size_t N, int S3, bool C3> 
    inline const SmallVectorView<T,M,S3,C3>& operator-=(
	const SmallVectorView<T,M,S3,C3>& v,
	const ProdMV<T,T1,T2>& pmv)
    {
      TMVAssert(pmv.size() == N);
      MultMV<true>(-pmv.GetX(), pmv.GetM(), pmv.GetV(), v); 
      return v; 
    }

  template <class T, size_t M, size_t N, int S3, bool C3> 
    inline const SmallVectorView<CT,M,S3,C3>& operator-=(
	const SmallVectorView<CT,M,S3,C3>& v,
	const ProdMV<T,T,T>& pmv)
    {
      TMVAssert(pmv.size() == N);
      MultMV<true>(-pmv.GetX(),pmv.GetM(),pmv.GetV(),v); 
      return v; 
    }

#define GENMATRIX GenSmallMatrix
#define GENVECTOR GenSmallVector
#define PRODMV Prodmv
#define PRODMV_1 Prodmv_1
#define PRODXM ProdXm
#define PRODXV ProdXv
#define X1a ,M,N,Si1,Sj1,C1
#define X1b ,N,M,Sj1,Si1,C1
#define X2 ,N,S2,C2
#define X3 ,M,N,Si1,Sj1,C1,S2,C2
#define Y ,size_t M,size_t N,int Si1,int Sj1,bool C1,int S2,bool C2
#include "TMV_AuxProdMV.h"
  // Defines things like v*m, m*v, (x*m)*(x*v), etc.
#define PRODMV_1 Prodmv_1
#define X3 ,M,N,Si1,Sj1,C1,S2,C2
#define Y ,size_t M,size_t N,int Si1,int Sj1,bool C1,int S2,bool C2
#include "TMV_AuxProdMVa.h"
  // Defines things like -(m*v), x*(m*v), (m*v)/x, etc.
#undef GENMATRIX
#undef GENVECTOR
#undef PRODMV
#undef PRODXM
#undef PRODXV

} // namespace tmv

#undef CT
#undef CCT
#undef VCT

#endif
