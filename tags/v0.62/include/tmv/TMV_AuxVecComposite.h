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


// This file sets up the Composite classes for all operations with a 
// (sparse) matrix that returns a vector
//
// Need to define the following with #define statements.
// (The given definition is for a Band Matrix.  Modify as 
// appropriate for the various other matrices.)
//
// #define GENMATRIX GenBandMatrix
// #define PRODMV ProdBV
// #define QUOTVM QuotVB
// #define RQUOTVM RQuotVB

//
// Matrix * Vector
//

template <class T, class T1, class T2> 
class PRODMV :
  public VectorComposite<T>
{
public:
  inline PRODMV(const T _x, const GENMATRIX<T1>& _m,
      const GenVector<T2>& _v) : x(_x), m(_m), v(_v)
  { TMVAssert(v.size()==m.rowsize()); }
  inline size_t size() const { return m.colsize(); }
  inline T GetX() const { return x; }
  inline const GENMATRIX<T1>& GetM() const { return m; }
  inline const GenVector<T2>& GetV() const { return v; }
  inline void AssignToV(const VectorView<RealType(T)>& v0) const
  {
    TMVAssert(IsReal(T()));
    TMVAssert(v0.size() == size());
    MultMV<false>(x,m,v,v0);
  }
  inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
  {
    TMVAssert(v0.size() == size());
    MultMV<false>(ComplexType(T)(x),m,v,v0);
  }
private:
  const T x;
  const GENMATRIX<T1>& m;
  const GenVector<T2>& v;
};

template <class T, class T1, class T2> 
inline const VectorView<T>& operator+=(
    const VectorView<T>& v, const PRODMV<T,T1,T2>& pmv)
{ MultMV<true>(pmv.GetX(),pmv.GetM(),pmv.GetV(),v); return v; }

template <class T> 
inline const VectorView<CT>& operator+=(
    const VectorView<CT>& v, const PRODMV<T,T,T>& pmv)
{ MultMV<true>(CT(pmv.GetX()),pmv.GetM(),pmv.GetV(),v); return v; }

template <class T, class T1, class T2> 
inline const VectorView<T>& operator-=(
    const VectorView<T>& v, const PRODMV<T,T1,T2>& pmv)
{ MultMV<true>(-pmv.GetX(), pmv.GetM(), pmv.GetV(), v); return v; }

template <class T> 
inline const VectorView<CT>& operator-=(
    const VectorView<CT>& v, const PRODMV<T,T,T>& pmv)
{ MultMV<true>(CT(-pmv.GetX()),pmv.GetM(),pmv.GetV(),v); return v; }


template <class T, class T1, class T2> 
class PRODVM :
  public VectorComposite<T>
{
public:
  inline PRODVM(const T _x, const GenVector<T1>& _v,
      const GENMATRIX<T2>& _m) : x(_x), v(_v), m(_m)
  { TMVAssert(v.size()==m.colsize()); }
  inline size_t size() const { return m.rowsize(); }
  inline T GetX() const { return x; }
  inline const GenVector<T1>& GetV() const { return v; }
  inline const GENMATRIX<T2>& GetM() const { return m; }
  inline void AssignToV(const VectorView<RealType(T)>& v0) const
  {
    TMVAssert(IsReal(T()));
    TMVAssert(v0.size() == size());
    MultMV<false>(x,m.Transpose(),v,v0);
  }
  inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
  {
    TMVAssert(v0.size() == size());
    MultMV<false>(ComplexType(T)(x),m.Transpose(),v,v0);
  }
private:
  const T x;
  const GenVector<T1>& v;
  const GENMATRIX<T2>& m;
};

template <class T> 
inline const VectorView<T>& operator*=(
    const VectorView<T>& v, const GENMATRIX<T>& m)
{ MultMV<false>(T(1),m.Transpose(),v,v); return v; }

template <class T> 
inline const VectorView<CT>& operator*=(
    const VectorView<CT>& v, const GENMATRIX<T>& m)
{ MultMV<false>(T(1),m.Transpose(),v,v); return v; }

template <class T, class T1, class T2> 
inline const VectorView<T>& operator+=(
    const VectorView<T>& v, const PRODVM<T,T1,T2>& pmv)
{ MultMV<true>(pmv.GetX(),pmv.GetM().Transpose(),pmv.GetV(),v); return v; }

template <class T> 
inline const VectorView<CT>& operator+=(
    const VectorView<CT>& v, const PRODVM<T,T,T>& pmv)
{ MultMV<true>(pmv.GetX(),pmv.GetM().Transpose(),pmv.GetV(),v); return v; }

template <class T, class T1, class T2> 
inline const VectorView<T>& operator-=(
    const VectorView<T>& v, const PRODVM<T,T1,T2>& pmv)
{ MultMV<true>(-pmv.GetX(), pmv.GetM().Transpose(), pmv.GetV(), v); return v; }

template <class T> 
inline const VectorView<CT>& operator-=(
    const VectorView<CT>& v, const PRODVM<T,T,T>& pmv)
{ MultMV<true>(-pmv.GetX(),pmv.GetM().Transpose(),pmv.GetV(),v); return v; }


#define GENMATRIX1 GENMATRIX
#define GENMATRIX2 GenVector
#define PRODMM PRODMV
#define PRODXM1 PRODXM
#define PRODXM2 ProdXV
#define GETM1 .GetM()
#define GETM2 .GetV()
#include "tmv/TMV_AuxProdMM.h"
#define GETM1 .GetM()
#define GETM2 .GetV()
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2
#define GENMATRIX1 GenVector
#define GENMATRIX2 GENMATRIX
#define PRODMM PRODVM
#define PRODXM1 ProdXV
#define PRODXM2 PRODXM
#define GETM1 .GetV()
#define GETM2 .GetM()
#include "tmv/TMV_AuxProdMM.h"
#define GETM1 .GetV()
#define GETM2 .GetM()
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2

//
// Vector / % Matrix 
// v/m is the solution (x) of mx = v
// ie. / is really division from the left: x = m^-1 v
// Use % if you want division from the right (v m^-1)
//

template <class T, class T1, class T2> 
class QUOTVM : 
  public VectorComposite<T>
{
public:
  inline QUOTVM(const T _x, const GenVector<T1>& _v,
      const GENMATRIX<T2>& _m) :
    x(_x), v(_v), m(_m)
  { TMVAssert(v.size()==m.colsize()); }
  inline size_t size() const { return m.rowsize(); }
  inline T GetX() const { return x; }
  inline const GenVector<T1>& GetV() const { return v; }
  inline const GENMATRIX<T2>& GetM() const { return m; }
  inline void AssignToV(const VectorView<RealType(T)>& v0) const
  {
    TMVAssert(IsReal(T()));
    TMVAssert(v0.size() == size());
    m.LDiv(v,v0);
    if (x != T(1)) v0 *= REAL(x);
  }
  inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
  {
    TMVAssert(v0.size() == size());
    m.LDiv(v,v0);
    if (x != T(1)) v0 *= x;
  }
private:
  const T x;
  const GenVector<T1>& v;
  const GENMATRIX<T2>& m;
};

template <class T, class T1, class T2> 
class RQUOTVM : 
  public VectorComposite<T>
{
public:
  inline RQUOTVM(const T _x, const GenVector<T1>& _v,
      const GENMATRIX<T2>& _m) :
    x(_x), v(_v), m(_m)
  { TMVAssert(v.size()==m.rowsize()); }
  inline size_t size() const { return m.colsize(); }
  inline T GetX() const { return x; }
  inline const GenVector<T1>& GetV() const { return v; }
  inline const GENMATRIX<T2>& GetM() const { return m; }
  inline void AssignToV(const VectorView<RealType(T)>& v0) const
  {
    TMVAssert(IsReal(T()));
    TMVAssert(v0.size() == size());
    m.RDiv(v,v0);
    if (x != T(1)) v0 *= REAL(x);
  }
  inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
  {
    TMVAssert(v0.size() == size());
    m.RDiv(v,v0);
    if (x != T(1)) v0 *= x;
  }
private:
  const T x;
  const GenVector<T1>& v;
  const GENMATRIX<T2>& m;
};

template <class T> 
inline const VectorView<T>& operator/=(
    const VectorView<T>& v, const GENMATRIX<T>& m)
{ 
  TMVAssert(m.IsSquare());
  TMVAssert(m.rowsize() == v.size());
  m.LDivEq(v); 
  return v; 
}

template <class T> 
inline const VectorView<CT>& operator/=(
    const VectorView<CT>& v, const GENMATRIX<T>& m)
{
  TMVAssert(m.IsSquare());
  TMVAssert(m.rowsize() == v.size());
  m.LDivEq(v); 
  return v; 
}

template <class T> 
inline const VectorView<T>& operator%=(
    const VectorView<T>& v, const GENMATRIX<T>& m)
{
  TMVAssert(m.IsSquare());
  TMVAssert(m.rowsize() == v.size());
  m.RDivEq(v); 
  return v; 
}

template <class T> 
inline const VectorView<CT>& operator%=(
    const VectorView<CT>& v, const GENMATRIX<T>& m)
{ 
  TMVAssert(m.IsSquare());
  TMVAssert(m.rowsize() == v.size());
  m.RDivEq(v); 
  return v; 
}

template <class T, class Tm> 
inline const VectorView<T>& operator*=(
    const VectorView<T>& v, const QUOTXM<T,Tm>& qxm)
{
  TMVAssert(qxm.GetM().IsSquare());
  TMVAssert(qxm.GetM().rowsize() == v.size());
  qxm.GetM().RDivEq(v); 
  v *= qxm.GetX(); 
  return v; 
}

template <class T> 
inline const VectorView<CT>& operator*=(
    const VectorView<CT>& v, const QUOTXM<T,T>& qxm)
{
  TMVAssert(qxm.GetM().IsSquare());
  TMVAssert(qxm.GetM().rowsize() == v.size());
  qxm.GetM().RDivEq(v); 
  v *= qxm.GetX(); 
  return v; 
}

#define GENMATRIX1 GenVector
#define GENMATRIX2 GENMATRIX
#define PRODXM1 ProdXV
#define PRODXM2 PRODXM
#define QUOTMM QUOTVM
#define RQUOTMM RQUOTVM
#define GETM1 .GetV()
#define GETM2 .GetM()
#include "tmv/TMV_AuxQuotMM.h"
#define GETM1 .GetV()
#define GETM2 .GetM()
#include "tmv/TMV_AuxQuotMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTMM
#undef RQUOTMM
