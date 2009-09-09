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


#ifndef TMV_SmallVectorArith_H
#define TMV_SmallVectorArith_H

#include "tmv/TMV_SmallVectorArithFunc.h"
#include "tmv/TMV_VectorArith.h"

#define CT std::complex<T>
#define CCT ConjRef<std::complex<T> >
#define VCT VarConjRef<std::complex<T> >

namespace tmv {

  template <class T, int N> 
  class SmallVectorComposite 
  {
  public:

    inline SmallVectorComposite() {}
    inline SmallVectorComposite(const SmallVectorComposite<T,N>&) {}
    virtual inline ~SmallVectorComposite() {}

    inline size_t size() const { return N; }
    virtual void AssignTov(SmallVector<RealType(T),N,CStyle>& v) const = 0;
    virtual void AssignTov(SmallVector<ComplexType(T),N,CStyle>& v) const = 0;
    virtual void AssignTov(SmallVector<RealType(T),N,FortranStyle>& v) const = 0;
    virtual void AssignTov(SmallVector<ComplexType(T),N,FortranStyle>& v) const = 0;
    virtual void AssignToV(const VectorView<RealType(T)>& v) const = 0;
    virtual void AssignToV(const VectorView<ComplexType(T)>& v) const = 0;
  };


  //
  // Vector * / Scalar
  //

  template <class T, class T1, int N, IndexStyle I> 
  class ProdXv : 
    public SmallVectorComposite<T,N> 
  {
  public:
    inline ProdXv(T _x, const SmallVector<T1,N,I>& _v) : x(_x), v(_v) {}
    inline T GetX() const { return x; }
    inline const SmallVector<T1,N,I>& GetV() const { return v; }
    inline void AssignTov(SmallVector<RealType(T),N,CStyle>& v0) const
    { TMVAssert(IsReal(T())); MultXV<N>(x,v.cptr(),v0.ptr()); }
    inline void AssignTov(SmallVector<ComplexType(T),N,CStyle>& v0) const
    { MultXV<N>(x,v.cptr(),v0.ptr()); }
    inline void AssignTov(SmallVector<RealType(T),N,FortranStyle>& v0) const
    { TMVAssert(IsReal(T())); MultXV<N>(x,v.cptr(),v0.ptr()); }
    inline void AssignTov(SmallVector<ComplexType(T),N,FortranStyle>& v0) const
    { MultXV<N>(x,v.cptr(),v0.ptr()); }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    { 
      TMVAssert(IsReal(T()));
      if (v0.step() == 1) 
        MultXV<N>(x,v.cptr(),v0.ptr());
      else
        MultXV(x,v.View(),v0); 
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    { 
      if (v0.step() == 1) {
        MultXV<N>(x,v.cptr(),v0.ptr());
        if (v0.isconj()) v0.ConjugateSelf();
      }
      else
        MultXV(x,v.View(),v0);
    }
  private:
    const T x;
    const SmallVector<T1,N,I>& v;
  };

  template <class T, int N, IndexStyle I> 
  inline SmallVector<T,N,I>& operator*=(SmallVector<T,N,I>& v1, T x2)
  { MultXV<N>(x2,v1.ptr()); return v1; }

  template <class T, int N, IndexStyle I> 
  inline SmallVector<CT,N,I>& operator*=(SmallVector<CT,N,I>& v1, T x2)
  { MultXV<N>(x2,v1.ptr()); return v1; }

  template <class T, int N, IndexStyle I> 
  inline SmallVector<CT,N,I>& operator*=(SmallVector<CT,N,I>& v1, CCT x2)
  { MultXV<N>(CT(x2),v1.ptr()); return v1; }

  template <class T, int N, IndexStyle I> 
  inline SmallVector<CT,N,I>& operator*=(SmallVector<CT,N,I>& v1, VCT x2)
  { MultXV<N>(CT(x2),v1.ptr()); return v1; }

  template <class T, int N, IndexStyle I> 
  inline SmallVector<T,N,I>& operator/=(SmallVector<T,N,I>& v1, T x2)
  { MultXV<N>(T(1)/x2,v1.ptr()); return v1; }

  template <class T, int N, IndexStyle I> 
  inline SmallVector<CT,N,I>& operator/=(SmallVector<CT,N,I>& v1, T x2)
  { MultXV<N>(T(1)/x2,v1.ptr()); return v1; }

  template <class T, int N, IndexStyle I> 
  inline SmallVector<CT,N,I>& operator/=(SmallVector<CT,N,I>& v1, CCT x2)
  { MultXV<N>(T(1)/CT(x2),v1.ptr()); return v1; }

  template <class T, int N, IndexStyle I> 
  inline SmallVector<CT,N,I>& operator/=(SmallVector<CT,N,I>& v1, VCT x2)
  { MultXV<N>(T(1)/CT(x2),v1.ptr()); return v1; }

#define GENMATRIX SmallVector
#define PRODXM ProdXv
#define X ,N,I
#define Y ,int N, IndexStyle I
#define GETM .GetV()
#include "tmv/TMV_AuxProdXM.h"
#undef GENMATRIX
#undef PRODXM

  //
  // Vector + Vector
  //

  template <class T, class T1, class T2, int N, IndexStyle I1, IndexStyle I2> 
  class Sumvv_1_1 :
    public SmallVectorComposite<T,N> 
  {
  public:
    inline Sumvv_1_1(
        T DEBUGPARAM(_x1), const SmallVector<T1,N,I1>& _v1, 
        T DEBUGPARAM(_x2), const SmallVector<T2,N,I2>& _v2) :
      v1(_v1),v2(_v2)
    { TMVAssert(_x1 == T(1) && _x2 == T(1)); }
    inline const SmallVector<T1,N,I1>& GetV1() const { return v1; }
    inline const SmallVector<T2,N,I2>& GetV2() const { return v2; }
    inline void AssignTov(SmallVector<RealType(T),N,CStyle>& v0) const
    { TMVAssert(IsReal(T())); AddVV_1_1<N>(v1.cptr(),v2.cptr(),v0.ptr()); }
    inline void AssignTov(SmallVector<ComplexType(T),N,CStyle>& v0) const
    { AddVV_1_1<N>(v1.cptr(),v2.cptr(),v0.ptr()); }
    inline void AssignTov(SmallVector<RealType(T),N,FortranStyle>& v0) const
    { TMVAssert(IsReal(T())); AddVV_1_1<N>(v1.cptr(),v2.cptr(),v0.ptr()); }
    inline void AssignTov(SmallVector<ComplexType(T),N,FortranStyle>& v0) const
    { AddVV_1_1<N>(v1.cptr(),v2.cptr(),v0.ptr()); }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    {
      TMVAssert(IsReal(T()));
      if (v0.step() == 1)
        AddVV_1_1<N>(v1.cptr(),v2.cptr(),v0.ptr());
      else
        AddVV(T(1),v1.View(),T(1),v2.View(),v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    {
      if (v0.step() == 1) {
        AddVV_1_1<N>(v1.cptr(),v2.cptr(),v0.ptr());
        if (v0.isconj()) v0.ConjugateSelf();
      }
      else
        AddVV(T(1),v1.View(),T(1),v2.View(),v0);
    }
  private:
    const SmallVector<T1,N,I1>& v1;
    const SmallVector<T2,N,I2>& v2;
  };

  template <class T, class T1, class T2, int N, IndexStyle I1, IndexStyle I2> 
  class Sumvv_1_m1 : 
    public SmallVectorComposite<T,N> 
  {
  public:
    inline Sumvv_1_m1(
        T DEBUGPARAM(_x1), const SmallVector<T1,N,I1>& _v1, 
        T DEBUGPARAM(_x2), const SmallVector<T2,N,I2>& _v2) :
      v1(_v1),v2(_v2)
    { TMVAssert(_x1 == T(1) && _x2 == T(-1)); }
    inline const SmallVector<T1,N,I1>& GetV1() const { return v1; }
    inline const SmallVector<T2,N,I2>& GetV2() const { return v2; }
    inline void AssignTov(SmallVector<RealType(T),N,CStyle>& v0) const
    { TMVAssert(IsReal(T())); AddVV_1_m1<N>(v1.cptr(),v2.cptr(),v0.ptr()); }
    inline void AssignTov(SmallVector<ComplexType(T),N,CStyle>& v0) const
    { AddVV_1_m1<N>(v1.cptr(),v2.cptr(),v0.ptr()); }
    inline void AssignTov(SmallVector<RealType(T),N,FortranStyle>& v0) const
    { TMVAssert(IsReal(T())); AddVV_1_m1<N>(v1.cptr(),v2.cptr(),v0.ptr()); }
    inline void AssignTov(SmallVector<ComplexType(T),N,FortranStyle>& v0) const
    { AddVV_1_m1<N>(v1.cptr(),v2.cptr(),v0.ptr()); }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    { 
      TMVAssert(IsReal(T()));
      if (v0.step() == 1)
        AddVV_1_m1<N>(v1.cptr(),v2.cptr(),v0.ptr());
      else
        AddVV(T(1),v1.View(),T(-1),v2.View(),v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    { 
      if (v0.step() == 1) {
        AddVV_1_m1<N>(v1.cptr(),v2.cptr(),v0.ptr());
        if (v0.isconj()) v0.ConjugateSelf();
      }
      else
        AddVV(T(1),v1.View(),T(-1),v2.View(),v0);
    }
  private:
    const SmallVector<T1,N,I1>& v1;
    const SmallVector<T2,N,I2>& v2;
  };

  template <class T, class T1, class T2, int N, IndexStyle I1, IndexStyle I2> 
  class Sumvv_1_x : 
    public SmallVectorComposite<T,N> 
  {
  public:
    inline Sumvv_1_x(
        T DEBUGPARAM(_x1), const SmallVector<T1,N,I1>& _v1, 
        T _x2, const SmallVector<T2,N,I2>& _v2) :
      v1(_v1),x2(_x2),v2(_v2)
    { TMVAssert(_x1 == T(1)); }
    inline const SmallVector<T1,N,I1>& GetV1() const { return v1; }
    inline T GetX2() const { return x2; }
    inline const SmallVector<T2,N,I2>& GetV2() const { return v2; }
    inline void AssignTov(SmallVector<RealType(T),N,CStyle>& v0) const
    { TMVAssert(IsReal(T())); AddVV_1_x<N>(v1.cptr(),x2,v2.cptr(),v0.ptr()); }
    inline void AssignTov(SmallVector<ComplexType(T),N,CStyle>& v0) const
    { AddVV_1_x<N>(v1.cptr(),x2,v2.cptr(),v0.ptr()); }
    inline void AssignTov(SmallVector<RealType(T),N,FortranStyle>& v0) const
    { TMVAssert(IsReal(T())); AddVV_1_x<N>(v1.cptr(),x2,v2.cptr(),v0.ptr()); }
    inline void AssignTov(SmallVector<ComplexType(T),N,FortranStyle>& v0) const
    { AddVV_1_x<N>(v1.cptr(),x2,v2.cptr(),v0.ptr()); }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    { 
      TMVAssert(IsReal(T()));
      if (v0.step() == 1)
        AddVV_1_x<N>(v1.cptr(),x2,v2.cptr(),v0.ptr());
      else
        AddVV(T(1),v1.View(),x2,v2.View(),v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    {
      if (v0.step() == 1) {
        AddVV_1_x<N>(v1.cptr(),x2,v2.cptr(),v0.ptr());
        if (v0.isconj()) v0.ConjugateSelf();
      }
      else
        AddVV(T(1),v1.View(),x2,v2.View(),v0);
    }
  private:
    const SmallVector<T1,N,I1>& v1;
    const T x2;
    const SmallVector<T2,N,I2>& v2;
  };

  template <class T, class T1, class T2, int N, IndexStyle I1, IndexStyle I2> 
  class Sumvv_x_1 : 
    public SmallVectorComposite<T,N> 
  {
  public:
    inline Sumvv_x_1(
        T _x1, const SmallVector<T1,N,I1>& _v1, 
        T DEBUGPARAM(_x2), const SmallVector<T2,N,I2>& _v2) :
      x1(_x1),v1(_v1),v2(_v2)
    { TMVAssert(_x2 == T(1)); }
    inline T GetX1() const { return x1; }
    inline const SmallVector<T1,N,I1>& GetV1() const { return v1; }
    inline const SmallVector<T2,N,I2>& GetV2() const { return v2; }
    inline void AssignTov(SmallVector<RealType(T),N,CStyle>& v0) const
    { TMVAssert(IsReal(T())); AddVV_x_1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr()); }
    inline void AssignTov(SmallVector<ComplexType(T),N,CStyle>& v0) const
    { AddVV_x_1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr()); }
    inline void AssignTov(SmallVector<RealType(T),N,FortranStyle>& v0) const
    { TMVAssert(IsReal(T())); AddVV_x_1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr()); }
    inline void AssignTov(SmallVector<ComplexType(T),N,FortranStyle>& v0) const
    { AddVV_x_1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr()); }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    { 
      TMVAssert(IsReal(T()));
      if (v0.step() == 1)
        AddVV_x_1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
      else
        AddVV(x1,v1.View(),T(1),v2.View(),v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    { 
      if (v0.step() == 1) {
        AddVV_x_1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
        if (v0.isconj()) v0.ConjugateSelf();
      }
      else
        AddVV(x1,v1.View(),T(1),v2.View(),v0);
    }
  private:
    const T x1;
    const SmallVector<T1,N,I1>& v1;
    const SmallVector<T2,N,I2>& v2;
  };

  template <class T, class T1, class T2, int N, IndexStyle I1, IndexStyle I2> 
  class Sumvv_x_m1 : 
    public SmallVectorComposite<T,N> 
  {
  public:
    inline Sumvv_x_m1(
        T _x1, const SmallVector<T1,N,I1>& _v1, 
        T DEBUGPARAM(_x2), const SmallVector<T2,N,I2>& _v2) :
      x1(_x1),v1(_v1),v2(_v2)
    { TMVAssert(_x2 == T(-1)); }
    inline T GetX1() const { return x1; }
    inline const SmallVector<T1,N,I1>& GetV1() const { return v1; }
    inline const SmallVector<T2,N,I2>& GetV2() const { return v2; }
    inline void AssignTov(
        SmallVector<RealType(T),N,CStyle>& v0) const
    { 
      TMVAssert(IsReal(T()));
      AddVV_x_m1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
    }
    inline void AssignTov(
        SmallVector<ComplexType(T),N,CStyle>& v0) const
    { 
      AddVV_x_m1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
    }
    inline void AssignTov(
        SmallVector<RealType(T),N,FortranStyle>& v0) const
    { 
      TMVAssert(IsReal(T()));
      AddVV_x_m1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
    }
    inline void AssignTov(
        SmallVector<ComplexType(T),N,FortranStyle>& v0) const
    {
      AddVV_x_m1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
    }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    {
      TMVAssert(IsReal(T()));
      if (v0.step() == 1)
        AddVV_x_m1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
      else
        AddVV(x1,v1.View(),T(-1),v2.View(),v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    {
      if (v0.step() == 1) {
        AddVV_x_m1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
        if (v0.isconj()) v0.ConjugateSelf();
      }
      else
        AddVV(x1,v1.View(),T(-1),v2.View(),v0);
    }
  private:
    const T x1;
    const SmallVector<T1,N,I1>& v1;
    const SmallVector<T2,N,I2>& v2;
  };

  template <class T, class T1, class T2, int N, IndexStyle I1, IndexStyle I2> 
  class Sumvv : 
    public SmallVectorComposite<T,N> 
  {
  public:
    inline Sumvv(T _x1, const SmallVector<T1,N,I1>& _v1, 
        T _x2, const SmallVector<T2,N,I2>& _v2) :
      x1(_x1),v1(_v1),x2(_x2), v2(_v2) {}
    inline T GetX1() const { return x1; }
    inline const SmallVector<T1,N,I1>& GetV1() const { return v1; }
    inline T GetX2() const { return x2; }
    inline const SmallVector<T2,N,I2>& GetV2() const { return v2; }
    inline void AssignTov(
        SmallVector<RealType(T),N,CStyle>& v0) const
    { 
      TMVAssert(IsReal(T()));
      AddVV<N>(x1,v1.cptr(),x2,v2.cptr(),v0.ptr());
    }
    inline void AssignTov(
        SmallVector<ComplexType(T),N,CStyle>& v0) const
    {
      AddVV<N>(x1,v1.cptr(),x2,v2.cptr(),v0.ptr());
    }
    inline void AssignTov(
        SmallVector<RealType(T),N,FortranStyle>& v0) const
    {
      TMVAssert(IsReal(T()));
      AddVV<N>(x1,v1.cptr(),x2,v2.cptr(),v0.ptr());
    }
    inline void AssignTov(
        SmallVector<ComplexType(T),N,FortranStyle>& v0) const
    {
      AddVV<N>(x1,v1.cptr(),x2,v2.cptr(),v0.ptr());
    }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    {
      TMVAssert(IsReal(T()));
      if (v0.step() == 1)
        AddVV<N>(x1,v1.cptr(),x2,v2.cptr(),v0.ptr());
      else
        AddVV(x1,v1.View(),x2,v2.View(),v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    {
      if (v0.step() == 1) {
        AddVV<N>(x1,v1.cptr(),x2,v2.cptr(),v0.ptr());
        if (v0.isconj()) v0.ConjugateSelf();
      }
      else
        AddVV(x1,v1.View(),x2,v2.View(),v0);
    }
  private:
    const T x1;
    const SmallVector<T1,N,I1>& v1;
    const T x2;
    const SmallVector<T2,N,I2>& v2;
  };

  template <class T, class T1, class T2, int N, IndexStyle I> 
  class SumVv : 
    public VectorComposite<T> 
  {
  public:
    inline SumVv(T _x1, const GenVector<T1>& _v1, 
        T _x2, const SmallVector<T2,N,I>& _v2) :
      x1(_x1),v1(_v1),x2(_x2), v2(_v2) { TMVAssert(v1.size() == N); }
    inline size_t size() const { return N; }
    inline T GetX1() const { return x1; }
    inline const GenVector<T1>& GetV1() const { return v1; }
    inline T GetX2() const { return x2; }
    inline const SmallVector<T2,N,I>& GetV2() const { return v2; }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    { 
      TMVAssert(IsReal(T()));
      if (v1.step() == 1 && v0.step() == 1 && (!SameStorage(v0,v1) )) {
        if (x1 == T(1)) 
          if (x2 == T(1))
            AddVV_1_1<N>(v1.cptr(),v2.cptr(),v0.ptr());
          else if (x2 == T(-1))
            AddVV_1_m1<N>(v1.cptr(),v2.cptr(),v0.ptr());
          else
            AddVV_1_x<N>(v1.cptr(),x2,v2.cptr(),v0.ptr());
        else 
          if (x2 == T(1))
            AddVV_x_1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
          else if (x2 == T(-1))
            AddVV_x_m1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
          else
            AddVV<N>(x1,v1.cptr(),x2,v2.cptr(),v0.ptr());
      }
      else AddVV(x1,v1,x2,v2.View(),v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    { 
      if (v1.step() == 1 && v0.step() == 1 && !v1.isconj() &&
          !SameStorage(v0,v1)) {
        if (x1 == T(1)) 
          if (x2 == T(1))
            AddVV_1_1<N>(v1.cptr(),v2.cptr(),v0.ptr());
          else if (x2 == T(-1))
            AddVV_1_m1<N>(v1.cptr(),v2.cptr(),v0.ptr());
          else
            AddVV_1_x<N>(v1.cptr(),x2,v2.cptr(),v0.ptr());
        else 
          if (x2 == T(1))
            AddVV_x_1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
          else if (x2 == T(-1))
            AddVV_x_m1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
          else
            AddVV<N>(x1,v1.cptr(),x2,v2.cptr(),v0.ptr());
        if (v0.isconj()) v0.ConjugateSelf();
      }
      else AddVV(x1,v1,x2,v2.View(),v0);
    }
  private:
    const T x1;
    const GenVector<T1>& v1;
    const T x2;
    const SmallVector<T2,N,I>& v2;
  };

  template <class T, class T1, class T2, int N, IndexStyle I> 
  class SumvV : 
    public VectorComposite<T> 
  {
  public:
    inline SumvV(T _x1, const SmallVector<T1,N,I>& _v1, 
        T _x2, const GenVector<T2>& _v2) :
      x1(_x1),v1(_v1),x2(_x2), v2(_v2) { TMVAssert(v2.size() == N); }
    inline size_t size() const { return N; }
    inline T GetX1() const { return x1; }
    inline const SmallVector<T1,N,I>& GetV1() const { return v1; }
    inline T GetX2() const { return x2; }
    inline const GenVector<T2>& GetV2() const { return v2; }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    { 
      TMVAssert(IsReal(T()));
      if (v2.step() == 1 && v0.step() == 1 && !SameStorage(v0,v2)) {
        if (x1 == T(1)) 
          if (x2 == T(1))
            AddVV_1_1<N>(v1.cptr(),v2.cptr(),v0.ptr());
          else if (x2 == T(-1))
            AddVV_1_m1<N>(v1.cptr(),v2.cptr(),v0.ptr());
          else
            AddVV_1_x<N>(v1.cptr(),x2,v2.cptr(),v0.ptr());
        else 
          if (x2 == T(1))
            AddVV_x_1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
          else if (x2 == T(-1))
            AddVV_x_m1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
          else
            AddVV<N>(x1,v1.cptr(),x2,v2.cptr(),v0.ptr());
      }
      else AddVV(x1,v1.View(),x2,v2,v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    { 
      if (v2.step() == 1 && v0.step() == 1 && !v2.isconj() &&
          !SameStorage(v0,v2) ) {
        if (x1 == T(1)) 
          if (x2 == T(1))
            AddVV_1_1<N>(v1.cptr(),v2.cptr(),v0.ptr());
          else if (x2 == T(-1))
            AddVV_1_m1<N>(v1.cptr(),v2.cptr(),v0.ptr());
          else
            AddVV_1_x<N>(v1.cptr(),x2,v2.cptr(),v0.ptr());
        else 
          if (x2 == T(1))
            AddVV_x_1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
          else if (x2 == T(-1))
            AddVV_x_m1<N>(x1,v1.cptr(),v2.cptr(),v0.ptr());
          else
            AddVV<N>(x1,v1.cptr(),x2,v2.cptr(),v0.ptr());
        if (v0.isconj()) v0.ConjugateSelf();
      }
      else AddVV(x1,v1.View(),x2,v2,v0);
    }
  private:
    const T x1;
    const SmallVector<T1,N,I>& v1;
    const T x2;
    const GenVector<T2>& v2;
  };

  // v+=v
  template <class T, int N, IndexStyle I1, IndexStyle I2> 
  inline SmallVector<T,N,I1>& operator+=(
      SmallVector<T,N,I1>& v1, const SmallVector<T,N,I2>& v2) 
  { AddVV_1<N>(v2.cptr(),v1.ptr()); return v1; }

  template <class T, int N, IndexStyle I1, IndexStyle I2> 
  inline SmallVector<CT,N,I1>& operator+=(
      SmallVector<CT,N,I1>& v1, const SmallVector<T,N,I2>& v2) 
  { AddVV_1<N>(v2.cptr(),v1.ptr()); return v1; }

  // v-=v
  template <class T, int N, IndexStyle I1, IndexStyle I2> 
  inline SmallVector<T,N,I1>& operator-=(
      SmallVector<T,N,I1>& v1, const SmallVector<T,N,I2>& v2)
  { AddVV_m1<N>(v2.cptr(),v1.ptr()); return v1; }

  template <class T, int N, IndexStyle I1, IndexStyle I2> 
  inline SmallVector<CT,N,I1>& operator-=(
      SmallVector<CT,N,I1>& v1, const SmallVector<T,N,I2>& v2) 
  { AddVV_m1<N>(v2.cptr(),v1.ptr()); return v1; }

  // v+=(x*v)
  template <class T, class T2, int N, IndexStyle I1, IndexStyle I2> 
  inline SmallVector<T,N,I1>& operator+=(
      SmallVector<T,N,I1>& v, const ProdXv<T,T2,N,I2>& v2)
  { AddVV<N>(v2.GetX(),v2.GetV().cptr(),v.ptr()); return v; }

  template <class T, int N, IndexStyle I1, IndexStyle I2> 
  inline SmallVector<CT,N,I1>& operator+=(
      SmallVector<CT,N,I1>& v, const ProdXv<T,T,N,I2>& v2)
  { AddVV<N>(v2.GetX(),v2.GetV().cptr(),v.ptr()); return v; }

  // v-=(x*v)
  template <class T, class T2, int N, IndexStyle I1, IndexStyle I2> 
  inline SmallVector<T,N,I1>& operator-=(
      SmallVector<T,N,I1>& v, const ProdXv<T,T2,N,I2>& v2)
  { AddVV<N>(-v2.GetX(),v2.GetV().cptr(),v.ptr()); return v; }

  template <class T, int N, IndexStyle I1, IndexStyle I2> 
  inline SmallVector<CT,N,I1>& operator-=(
      SmallVector<CT,N,I1>& v, const ProdXv<T,T,N,I2>& v2)
  { AddVV<N>(-v2.GetX(),v2.GetV().cptr(),v.ptr()); return v; }

  // Mix with Vector
  template <class T, int N, IndexStyle I> 
  inline SmallVector<T,N,I>& operator+=(
      SmallVector<T,N,I>& v1, const GenVector<T>& v2) 
  {
    //std::cout<<"v += V\n";
    //std::cout<<"v1 = "<<v1<<std::endl;
    //std::cout<<"v2 = "<<TypeText(v2)<<"  "<<v2<<std::endl;
    //std::cout<<"v1+v2 = "<<(tmv::Vector<T>(v1)+tmv::Vector<T>(v2))<<std::endl;
    TMVAssert(v2.size() == N);
    if (v2.step() == 1 && !v2.isconj())
      AddVV_1<N>(v2.cptr(),v1.ptr());
    else
      AddVV(T(1),v2,v1.View()); 
    //std::cout<<"v1 => "<<v1<<std::endl;
    return v1; 
  }

  template <class T, int N, IndexStyle I> 
  inline SmallVector<CT,N,I>& operator+=(
      SmallVector<CT,N,I>& v1, const GenVector<T>& v2) 
  {
    //std::cout<<"v += V\n";
    //std::cout<<"v1 = "<<v1<<std::endl;
    //std::cout<<"v2 = "<<TypeText(v2)<<"  "<<v2<<std::endl;
    //std::cout<<"v1+v2 = "<<(tmv::Vector<T>(v1)+tmv::Vector<T>(v2))<<std::endl;
    TMVAssert(v2.size() == N);
    if (v2.step() == 1)
      AddVV_1<N>(v2.cptr(),v1.ptr());
    else
      AddVV(T(1),v2,v1.View()); 
    return v1; 
    //std::cout<<"v1 => "<<v1<<std::endl;
  }

  template <class T, int N, IndexStyle I> 
  inline SmallVector<T,N,I>& operator-=(
      SmallVector<T,N,I>& v1, const GenVector<T>& v2) 
  {
    TMVAssert(v2.size() == N);
    if (v2.step() == 1 && !v2.isconj())
      AddVV_m1<N>(v2.cptr(),v1.ptr());
    else
      AddVV(T(-1),v2,v1.View()); 
    return v1; 
  }

  template <class T, int N, IndexStyle I> 
  inline SmallVector<CT,N,I>& operator-=(
      SmallVector<CT,N,I>& v1, const GenVector<T>& v2) 
  {
    TMVAssert(v2.size() == N);
    if (v2.step() == 1)
      AddVV_m1<N>(v2.cptr(),v1.ptr());
    else
      AddVV(T(-1),v2,v1.View()); 
    return v1; 
  }

  template <class T, class Tv> 
  class ProdXV;

  template <class T, class T2, int N, IndexStyle I> 
  inline SmallVector<T,N,I>& operator+=(
      SmallVector<T,N,I>& v, const ProdXV<T,T2>& v2) 
  {
    //std::cout<<"v += XV\n";
    //std::cout<<"v = "<<v<<std::endl;
    //std::cout<<"v2 = "<<TypeText(v2.GetV())<<"  "<<v2.GetV()<<std::endl;
    TMVAssert(v2.size() == N);
    if (v2.GetV().step() == 1 && !v2.GetV().isconj())
      AddVV<N>(v2.GetX(),v2.GetV().cptr(),v.ptr());
    else
      AddVV(v2.GetX(),v2.GetV(),v.View()); 
    //std::cout<<"v => "<<v1<<std::endl;
    return v; 
  }

  template <class T, int N, IndexStyle I> 
  inline SmallVector<CT,N,I>& operator+=(
      SmallVector<CT,N,I>& v, const ProdXV<T,T>& v2) 
  {
    //std::cout<<"v += XV\n";
    //std::cout<<"v = "<<v<<std::endl;
    //std::cout<<"v2 = "<<TypeText(v2.GetV())<<"  "<<v2.GetV()<<std::endl;
    TMVAssert(v2.size() == N);
    if (v2.GetV().step() == 1 && !v2.GetV().isconj())
      AddVV<N>(v2.GetX(),v2.GetV().cptr(),v.ptr());
    else
      AddVV(v2.GetX(),v2.GetV(),v.View()); 
    //std::cout<<"v => "<<v<<std::endl;
    return v; 
  }

  template <class T, class T2, int N, IndexStyle I> 
  inline SmallVector<T,N,I>& operator-=(
      SmallVector<T,N,I>& v, const ProdXV<T,T2>& v2) 
  {
    TMVAssert(v2.size() == N);
    if (v2.GetV().step() == 1 && !v2.GetV().isconj())
      AddVV<N>(-v2.GetX(),v2.GetV().cptr(),v.ptr());
    else
      AddVV(-v2.GetX(),v2.GetV(),v.View()); 
    return v; 
  }

  template <class T, int N, IndexStyle I> 
  inline SmallVector<CT,N,I>& operator-=(
      SmallVector<CT,N,I>& v, const ProdXV<T,T>& v2) 
  {
    TMVAssert(v2.size() == N);
    if (v2.GetV().step() == 1 && !v2.GetV().isconj())
      AddVV<N>(-v2.GetX(),v2.GetV().cptr(),v.ptr());
    else
      AddVV(-v2.GetX(),v2.GetV(),v.View()); 
    return v; 
  }

  template <class T, int N, IndexStyle I> 
  inline const VectorView<T>& operator+=(
      const VectorView<T>& v1, const SmallVector<T,N,I>& v2) 
  { 
    //std::cout<<"V += v\n";
    //std::cout<<"v1 = "<<TypeText(v1)<<"  "<<v1<<std::endl;
    //std::cout<<"v2 = "<<v2<<std::endl;
    if (v1.step() == 1 && !v1.isconj())
      AddVV_1<N>(v2.cptr(),v1.ptr());
    else
      AddVV(T(1),v2.View(),v1);
    //std::cout<<"v1 => "<<v1<<std::endl;
    return v1; 
  }

  template <class T, int N, IndexStyle I> 
  inline const VectorView<CT>& operator+=(
      const VectorView<CT>& v1, const SmallVector<T,N,I>& v2) 
  {
    //std::cout<<"V += v\n";
    //std::cout<<"v1 = "<<TypeText(v1)<<"  "<<v1<<std::endl;
    //std::cout<<"v2 = "<<v2<<std::endl;
    if (v1.step() == 1)
      AddVV_1<N>(v2.cptr(),v1.ptr());
    else
      AddVV(T(1),v2.View(),v1);
    //std::cout<<"v1 => "<<v1<<std::endl;
    return v1;
  }

  // v-=v
  template <class T, int N, IndexStyle I> 
  inline const VectorView<T>& operator-=(
      const VectorView<T>& v1, const SmallVector<T,N,I>& v2)
  { 
    if (v1.step() == 1 && !v1.isconj())
      AddVV_m1<N>(v2.cptr(),v1.ptr());
    else
      AddVV(T(-1),v2.View(),v1);
    return v1;
  }

  template <class T, int N, IndexStyle I> 
  inline const VectorView<CT>& operator-=(
      const VectorView<CT>& v1, const SmallVector<T,N,I>& v2) 
  { 
    if (v1.step() == 1)
      AddVV_m1<N>(v2.cptr(),v1.ptr());
    else
      AddVV(T(-1),v2.View(),v1);
    return v1;
  }

  // v+=(x*v)
  template <class T, class T2, int N, IndexStyle I> 
  inline const VectorView<T>& operator+=(
      const VectorView<T>& v, const ProdXv<T,T2,N,I>& v2)
  {
    if (v.step() == 1 && !v.isconj())
      AddVV<N>(v2.GetX(),v2.GetV().cptr(),v.ptr());
    else
      AddVV(v2.GetX(),v2.GetV().View(),v);
    return v;
  }

  template <class T, int N, IndexStyle I> 
  inline const VectorView<CT>& operator+=(
      const VectorView<CT>& v, const ProdXv<T,T,N,I>& v2)
  { 
    if (v.step() == 1)
      AddVV<N>(v2.GetX(),v2.GetV().cptr(),v.ptr());
    else
      AddVV(v2.GetX(),v2.GetV().View(),v);
    return v;
  }

  // v-=(x*v)
  template <class T, class T2, int N, IndexStyle I> 
  inline const VectorView<T>& operator-=(
      const VectorView<T>& v, const ProdXv<T,T2,N,I>& v2)
  { 
    if (v.step() == 1 && !v.isconj())
      AddVV<N>(-v2.GetX(),v2.GetV().cptr(),v.ptr());
    else
      AddVV(-v2.GetX(),v2.GetV().View(),v);
    return v;
  }

  template <class T, int N, IndexStyle I> 
  inline const VectorView<CT>& operator-=(
      const VectorView<CT>& v, const ProdXv<T,T,N,I>& v2)
  {
    if (v.step() == 1)
      AddVV<N>(-v2.GetX(),v2.GetV().cptr(),v.ptr());
    else
      AddVV(-v2.GetX(),v2.GetV().View(),v);
    return v;
  }

#define GENMATRIX1 SmallVector
#define GENMATRIX2 SmallVector
#define SUMMM Sumvv
#define SUMMM_1_1 Sumvv_1_1
#define SUMMM_1_m1 Sumvv_1_m1
#define SUMMM_1_x Sumvv_1_x
#define SUMMM_x_1 Sumvv_x_1
#define SUMMM_x_m1 Sumvv_x_m1
#define PRODXM1 ProdXv
#define PRODXM2 ProdXv
#define X1 ,N,I1
#define X2 ,N,I2
#define X3 ,N,I1,I2
#define Y ,int N, IndexStyle I1, IndexStyle I2
#define GETM .GetV()
#include "tmv/TMV_AuxSumMM.h"
  // These get undef'ed in TMV_AuxSumMM.h
#define SUMMM_1_1 Sumvv_1_1
#define SUMMM_1_m1 Sumvv_1_m1
#define SUMMM_1_x Sumvv_1_x
#define SUMMM_x_1 Sumvv_x_1
#define SUMMM_x_m1 Sumvv_x_m1
#define X3 ,N,I1,I2
#define Y ,int N, IndexStyle I1, IndexStyle I2
#define GETM .GetV()
#include "tmv/TMV_AuxSumMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef SUMMM
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 SmallVector
#define GENMATRIX2 GenVector
#define SUMMM SumvV
#define PRODXM1 ProdXv
#define PRODXM2 ProdXV
#define X1 ,N,I
#define X3 ,N,I
#define Y ,int N, IndexStyle I
#define GETM .GetV()
#include "tmv/TMV_AuxSumMM.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef SUMMM
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 GenVector
#define GENMATRIX2 SmallVector
#define SUMMM SumVv
#define PRODXM1 ProdXV
#define PRODXM2 ProdXv
#define X2 ,N,I
#define X3 ,N,I
#define Y ,int N, IndexStyle I
#define GETM .GetV()
#include "tmv/TMV_AuxSumMM.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef SUMMM
#undef PRODXM1
#undef PRODXM2

  //
  // Vector * Vector
  //

  template <class T, int N, IndexStyle I1, IndexStyle I2> 
  inline T operator*(
      const SmallVector<T,N,I1>& v1, const SmallVector<T,N,I2>& v2) 
  { return MultVV<N>(v1.cptr(),v2.cptr()); }

  template <class T, int N, IndexStyle I1, IndexStyle I2> 
  inline CT operator*(
      const SmallVector<CT,N,I1>& v1, const SmallVector<T,N,I2>& v2)
  { return MultVV<N>(v1.cptr(),v2.cptr()); }

  template <class T, int N, IndexStyle I1, IndexStyle I2> 
  inline CT operator*(
      const SmallVector<T,N,I1>& v1, const SmallVector<CT,N,I2>& v2)
  { return MultVV<N>(v2.cptr(),v1.cptr()); }

  // v * (x*v)
  template <class T, class T2, int N, IndexStyle I1, IndexStyle I2> 
  inline T operator*(
      const SmallVector<T,N,I1>& v1, const ProdXv<T,T2,N,I2>& v2) 
  { return v2.GetX()*MultVV<N>(v1.cptr(),v2.GetV().cptr()); }

  template <class T, int N, IndexStyle I1, IndexStyle I2> 
  inline CT operator*(
      const SmallVector<CT,N,I1>& v1, const ProdXv<T,T,N,I2>& v2)
  { return v2.GetX()*MultVV<N>(v1.cptr(),v2.GetV().cptr()); }

  template <class T, class T2, int N, IndexStyle I1, IndexStyle I2> 
  inline CT operator*(
      const SmallVector<T,N,I1>& v1, const ProdXv<CT,T2,N,I2>& v2)
  { return v2.GetX()*MultVV<N>(v2.GetV().cptr(),v1.cptr()); }

  // (x*v) * v
  template <class T, class T1, int N, IndexStyle I1, IndexStyle I2> 
  inline T operator*(
      const ProdXv<T,T1,N,I1>& v1, const SmallVector<T,N,I2>& v2)
  { return v1.GetX()*MultVV<N>(v1.GetV().cptr(),v2.cptr()); }

  template <class T, int N, IndexStyle I1, IndexStyle I2> 
  inline CT operator*(
      const ProdXv<T,T,N,I1>& v1, const SmallVector<CT,N,I2>& v2)
  { return v1.GetX()*MultVV<N>(v2.cptr(),v1.GetV().cptr()); }

  template <class T, class T1, int N, IndexStyle I1, IndexStyle I2> 
  inline CT operator*(
      const ProdXv<CT,T1,N,I1>& v1, const SmallVector<T,N,I2>& v2)
  { return v1.GetX()*MultVV<N>(v1.GetV().cptr(),v2.cptr()); }

  // (x*v) * (x*v)
  template <class T, class T1, class T2, int N, IndexStyle I1, IndexStyle I2> 
  inline T operator*(
      const ProdXv<T,T1,N,I1>& v1, const ProdXv<T,T2,N,I2>& v2)
  { return v1.GetX()*v2.GetX()*MultVV<N>(v1.GetV().cptr(),v2.GetV().cptr()); }

  template <class T, class T1, int N, IndexStyle I1, IndexStyle I2> 
  inline CT operator*(
      const ProdXv<CT,T1,N,I1>& v1, const ProdXv<T,T,N,I2>& v2)
  { return v1.GetX()*v2.GetX()*MultVV<N>(v1.GetV().cptr(),v2.GetV().cptr()); }

  template <class T, class T2, int N, IndexStyle I1, IndexStyle I2> 
  inline CT operator*(
      const ProdXv<T,T,N,I1>& v1, const ProdXv<CT,T2,N,I2>& v2)
  { return v1.GetX()*v2.GetX()*MultVV<N>(v1.GetV().cptr(),v2.GetV().cptr()); }


  // Mix with Vector:

  // v * v
  template <class T, int N, IndexStyle I> 
  inline T operator*(const GenVector<T>& v1, const SmallVector<T,N,I>& v2) 
  { 
    if (v1.step() == 1 && !v1.isconj())
      return MultVV<N>(v1.cptr(),v2.cptr());
    else 
      return MultVV(v1,v2.View());
  }

  template <class T, int N, IndexStyle I> 
  inline CT operator*(const GenVector<CT>& v1, const SmallVector<T,N,I>& v2)
  { 
    if (v1.step() == 1 && !v1.isconj())
      return MultVV<N>(v1.cptr(),v2.cptr());
    else 
      return MultVV(v1,v2.View()); 
  }

  template <class T, int N, IndexStyle I> 
  inline CT operator*(const GenVector<T>& v1, const SmallVector<CT,N,I>& v2)
  {
    if (v1.step() == 1 && !v1.isconj())
      return MultVV<N>(v2.cptr(),v1.cptr());
    else 
      return MultVV(v2.View(),v1);
  }

  // v * (x*v)
  template <class T, class T2, int N, IndexStyle I> 
  inline T operator*(const GenVector<T>& v1, const ProdXv<T,T2,N,I>& v2) 
  { 
    if (v1.step() == 1 && !v1.isconj())
      return v2.GetX()*MultVV<N>(v1.cptr(),v2.GetV().cptr());
    else 
      return v2.GetX()*MultVV(v1,v2.GetV().View());
  }

  template <class T, int N, IndexStyle I> 
  inline CT operator*(const GenVector<CT>& v1, const ProdXv<T,T,N,I>& v2)
  {
    if (v1.step() == 1 && !v1.isconj())
      return v2.GetX()*MultVV<N>(v1.cptr(),v2.GetV().cptr());
    else 
      return v2.GetX()*MultVV(v1,v2.GetV().View());
  }

  template <class T, class T2, int N, IndexStyle I> 
  inline CT operator*(const GenVector<T>& v1, const ProdXv<CT,T2,N,I>& v2)
  {
    if (v1.step() == 1 && !v1.isconj())
      return v2.GetX()*MultVV<N>(v2.GetV().cptr(),v1.cptr());
    else 
      return v2.GetX()*MultVV(v2.GetV().View(),v1);
  }

  // (x*v) * v
  template <class T, class T1, int N, IndexStyle I> 
  inline T operator*(const ProdXV<T,T1>& v1, const SmallVector<T,N,I>& v2)
  {
    if (v1.GetV().step() == 1)
      return v1.GetX()*MultVV<N>(v1.GetV().cptr(),v2.cptr());
    else 
      return v1.GetX()*MultVV(v1.GetV(),v2.View());
  }

  template <class T, class T1, int N, IndexStyle I> 
  inline CT operator*(const ProdXV<CT,T1>& v1, const SmallVector<T,N,I>& v2)
  {
    if (v1.GetV().step() == 1)
      return v1.GetX()*MultVV<N>(v1.GetV().cptr(),v2.cptr());
    else 
      return v1.GetX()*MultVV(v1.GetV(),v2.View());
  }

  template <class T, int N, IndexStyle I> 
  inline CT operator*(const ProdXV<T,T>& v1, const SmallVector<CT,N,I>& v2)
  {
    if (v1.GetV().step() == 1)
      return v1.GetX()*MultVV<N>(v2.cptr(),v1.GetV().cptr());
    else 
      return v1.GetX()*MultVV(v2.View(),v1.GetV());
  }

  // (x*v) * (x*v)
  template <class T, class T1, class T2, int N, IndexStyle I> 
  inline T operator*(const ProdXV<T,T1>& v1, const ProdXv<T,T2,N,I>& v2)
  {
    if (v1.GetV().step() == 1)
      return v1.GetX()*v2.GetX()*MultVV<N>(v1.GetV().cptr(),v2.GetV().cptr());
    else 
      return v1.GetX()*v2.GetX()*MultVV(v1.GetV(),v2.GetV().View());
  }

  template <class T, class T1, int N, IndexStyle I> 
  inline CT operator*(const ProdXV<CT,T1>& v1, const ProdXv<T,T,N,I>& v2)
  {
    if (v1.GetV().step() == 1)
      return v1.GetX()*v2.GetX()*MultVV<N>(v1.GetV().cptr(),v2.GetV().cptr());
    else 
      return v1.GetX()*v2.GetX()*MultVV(v1.GetV(),v2.GetV().View());
  }

  template <class T, class T2, int N, IndexStyle I> 
  inline CT operator*(const ProdXV<T,T>& v1, const ProdXv<CT,T2,N,I>& v2)
  {
    if (v1.GetV().step() == 1)
      return v1.GetX()*v2.GetX()*MultVV<N>(v2.GetV().cptr(),v1.GetV().cptr());
    else 
      return v1.GetX()*v2.GetX()*MultVV(v2.GetV().View(),v1.GetV());
  }

  // v * v 
  template <class T, int N, IndexStyle I> 
  inline T operator*(const SmallVector<T,N,I>& v1, const GenVector<T>& v2) 
  {
    if (v2.step() == 1 && !v2.step())
      return MultVV<N>(v1.cptr(),v2.cptr());
    else 
      return MultVV(v1.View(),v2);
  }

  template <class T, int N, IndexStyle I> 
  inline CT operator*(const SmallVector<CT,N,I>& v1, const GenVector<T>& v2)
  {
    if (v2.step() == 1 && !v2.step())
      return MultVV<N>(v1.cptr(),v2.cptr());
    else 
      return MultVV(v1.View(),v2);
  }

  template <class T, int N, IndexStyle I> 
  inline CT operator*(const SmallVector<T,N,I>& v1, const GenVector<CT>& v2)
  { 
    if (v2.step() == 1 && !v2.step())
      return MultVV<N>(v2.cptr(),v1.cptr());
    else 
      return MultVV(v2,v1.View());
  }

  // v * (x*v)
  template <class T, class T2, int N, IndexStyle I> 
  inline T operator*(const SmallVector<T,N,I>& v1, const ProdXV<T,T2>& v2) 
  { 
    if (v2.GetV().step() == 1)
      return v2.GetX()*MultVV<N>(v1.cptr(),v2.GetV().cptr());
    else 
      return v2.GetX()*MultVV(v1.View(),v2.GetV());
  }

  template <class T, int N, IndexStyle I> 
  inline CT operator*(const SmallVector<CT,N,I>& v1, const ProdXV<T,T>& v2)
  {
    if (v2.GetV().step() == 1)
      return v2.GetX()*MultVV<N>(v1.cptr(),v2.GetV().cptr());
    else 
      return v2.GetX()*MultVV(v1.View(),v2.GetV());
  }

  template <class T, class T2, int N, IndexStyle I> 
  inline CT operator*(const SmallVector<T,N,I>& v1, const ProdXV<CT,T2>& v2)
  { 
    if (v2.GetV().step() == 1)
      return v2.GetX()*MultVV<N>(v2.GetV().cptr(),v1.cptr());
    else 
      return v2.GetX()*MultVV(v2.GetV(),v1.View());
  }

  // (x*v) * v
  template <class T, class T1, int N, IndexStyle I> 
  inline T operator*(const ProdXv<T,T1,N,I>& v1, const GenVector<T>& v2)
  {
    if (v2.step() == 1 && !v2.step())
      return v1.GetX()*MultVV<N>(v1.GetV().cptr(),v2.cptr());
    else 
      return v1.GetX()*MultVV(v1.GetV().View(),v2);
  }

  template <class T, class T1, int N, IndexStyle I> 
  inline CT operator*(const ProdXv<CT,T1,N,I>& v1, const GenVector<T>& v2)
  {
    if (v2.step() == 1 && !v2.step())
      return v1.GetX()*MultVV<N>(v1.GetV().cptr(),v2.cptr());
    else 
      return v1.GetX()*MultVV(v1.GetV().View(),v2);
  }

  template <class T, int N, IndexStyle I> 
  inline CT operator*(const ProdXv<T,T,N,I>& v1, const GenVector<CT>& v2)
  {
    if (v2.step() == 1 && !v2.step())
      return v1.GetX()*MultVV<N>(v2.cptr(),v1.GetV().cptr());
    else 
      return v1.GetX()*MultVV(v2,v1.GetV().View());
  }

  // (x*v) * (x*v)

  template <class T, class T1, class T2, int N, IndexStyle I> 
  inline T operator*(const ProdXv<T,T1,N,I>& v1, const ProdXV<T,T2>& v2)
  {
    if (v2.GetV().step() == 1)
      return v1.GetX()*v2.GetX()*MultVV<N>(v1.GetV().cptr(),v2.GetV().cptr());
    else 
      return v1.GetX()*v2.GetX()*MultVV(v1.GetV().View(),v2.GetV());
  }

  template <class T, class T1, int N, IndexStyle I> 
  inline CT operator*(const ProdXv<CT,T1,N,I>& v1, const ProdXV<T,T>& v2)
  {
    if (v2.GetV().step() == 1)
      return v1.GetX()*v2.GetX()*MultVV<N>(v1.GetV().cptr(),v2.GetV().cptr());
    else 
      return v1.GetX()*v2.GetX()*MultVV(v1.GetV().View(),v2.GetV());
  }

  template <class T, class T2, int N, IndexStyle I> 
  inline CT operator*(const ProdXv<T,T,N,I>& v1, const ProdXV<CT,T2>& v2)
  {
    if (v2.GetV().step() == 1)
      return v1.GetX()*v2.GetX()*MultVV<N>(v2.GetV().cptr(),v1.GetV().cptr());
    else 
      return v1.GetX()*v2.GetX()*MultVV(v2.GetV(),v1.GetV().View());
  }


} // namespace tmv

#undef CT
#undef CCT
#undef VCT

#endif 
