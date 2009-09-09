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


#ifndef TMV_SmallMatrixArith_H
#define TMV_SmallMatrixArith_H

#include "tmv/TMV_SmallVectorArithFunc.h"
#include "tmv/TMV_SmallMatrixArithFunc.h"
#include "tmv/TMV_MatrixArith.h"

#define CT std::complex<T>
#define CCT ConjRef<std::complex<T> >
#define VCT VarConjRef<std::complex<T> >

namespace tmv {

  template <class T, int M, int N> 
  class SmallMatrixComposite 
  {
  public:

    inline SmallMatrixComposite() {}
    inline SmallMatrixComposite(const SmallMatrixComposite<T,M,N>&) {}
    virtual inline ~SmallMatrixComposite() {}

    size_t colsize() const { return M; }
    size_t rowsize() const { return N; }

    virtual void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,CStyle>& m) const = 0;
    virtual void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,CStyle>& m) const = 0;
    virtual void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,CStyle>& m) const = 0;
    virtual void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,CStyle>& m) const = 0;
    virtual void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,FortranStyle>& m) const = 0;
    virtual void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,FortranStyle>& m) const = 0;
    virtual void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,FortranStyle>& m) const = 0;
    virtual void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,FortranStyle>& m) const = 0;
    virtual void AssignToM(const MatrixView<RealType(T)>& m) const = 0;
    virtual void AssignToM(const MatrixView<ComplexType(T)>& m) const = 0;
  };


  //
  // Scalar * Matrix
  //

  template <class T, class T1, int M, int N, StorageType S, IndexStyle I> 
  class ProdXm : 
    public SmallMatrixComposite<T,M,N> 
  {
  public:
    inline ProdXm(const T _x, const SmallMatrix<T1,M,N,S,I>& _m) :
      x(_x), m(_m) {}
    inline T GetX() const { return x; }
    inline const SmallMatrix<T1,M,N,S,I>& GetM() const { return m; }
    inline StorageType stor() const { return S; }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,CStyle>& m0) const
    { 
      TMVAssert(IsReal(T()));
      m0=m;
      MultXV<M*N>(x,m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,CStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0=m;
      MultXV<M*N>(x,m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,CStyle>& m0) const
    {
      m0=m;
      MultXV<M*N>(x,m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,CStyle>& m0) const
    {
      m0=m;
      MultXV<M*N>(x,m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,FortranStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0=m;
      MultXV<M*N>(x,m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,FortranStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0=m;
      MultXV<M*N>(x,m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,FortranStyle>& m0) const
    {
      m0=m;
      MultXV<M*N>(x,m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,FortranStyle>& m0) const
    {
      m0=m;
      MultXV<M*N>(x,m0.ptr());
    }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    {
      TMVAssert(IsReal(T()));
      if (m0.stepi() == 1 && m0.stepj() == M && !SameStorage(m0,m)) {
        DoCopy<M,N,S,ColMajor>(m.cptr(),m0.ptr());
        MultXV<M*N>(x,m0.ptr());
      } else if (m0.stepj() == 1 && m0.stepi() == N && !SameStorage(m0,m)) {
        DoCopy<M,N,S,RowMajor>(m.cptr(),m0.ptr());
        MultXV<M*N>(x,m0.ptr());
      } else {
        MultXM(x,m0=m);
      }
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    {
      if (m0.stepi() == 1 && m0.stepj() == M && !SameStorage(m0,m)) {
        DoCopy<M,N,S,ColMajor>(m.cptr(),m0.ptr());
        MultXV<M*N>(x,m0.ptr());
        if (m0.isconj()) m0.ConjugateSelf();
      } else if (m0.stepj() == 1 && m0.stepi() == N && !SameStorage(m0,m)) {
        DoCopy<M,N,S,RowMajor>(m.cptr(),m0.ptr());
        MultXV<M*N>(x,m0.ptr());
        if (m0.isconj()) m0.ConjugateSelf();
      } else {
        MultXM(x,m0=m);
      }
    }
  private:
    const T x;
    const SmallMatrix<T1,M,N,S,I>& m;
  };

  // m*=x
  template <class T, int M, int N, StorageType S, IndexStyle I>
  inline SmallMatrix<T,M,N,S,I>& operator*=(SmallMatrix<T,M,N,S,I>& m, T x) 
  { MultXV<M*N>(x,m.ptr()); return m; }

  template <class T, int M, int N, StorageType S, IndexStyle I>
  inline SmallMatrix<CT,M,N,S,I>& operator*=(SmallMatrix<CT,M,N,S,I>& m, T x) 
  { MultXV<M*N>(x,m.ptr()); return m; }

  template <class T, int M, int N, StorageType S, IndexStyle I>
  inline SmallMatrix<CT,M,N,S,I>& operator*=(SmallMatrix<CT,M,N,S,I>& m, CCT x) 
  { m *= CT(x); return m; }

  template <class T, int M, int N, StorageType S, IndexStyle I>
  inline SmallMatrix<CT,M,N,S,I>& operator*=(SmallMatrix<CT,M,N,S,I>& m, VCT x) 
  { m *= CT(x); return m; }

  // m/=x
  template <class T, int M, int N, StorageType S, IndexStyle I>
  inline SmallMatrix<T,M,N,S,I>& operator/=(SmallMatrix<T,M,N,S,I>& m, T x) 
  { m *= T(1)/x; return m; }

  template <class T, int M, int N, StorageType S, IndexStyle I>
  inline SmallMatrix<CT,M,N,S,I>& operator/=(SmallMatrix<CT,M,N,S,I>& m, T x) 
  { m *= T(1)/x; return m; }

  template <class T, int M, int N, StorageType S, IndexStyle I>
  inline SmallMatrix<CT,M,N,S,I>& operator/=(SmallMatrix<CT,M,N,S,I>& m, CCT x) 
  { m *= T(1)/CT(x); return m; }

  template <class T, int M, int N, StorageType S, IndexStyle I>
  inline SmallMatrix<CT,M,N,S,I>& operator/=(SmallMatrix<CT,M,N,S,I>& m, VCT x) 
  { m *= T(1)/CT(x); return m; }

#define GENMATRIX SmallMatrix
#define PRODXM ProdXm
#define X ,M,N,S,I
#define Y ,int M, int N, StorageType S, IndexStyle I
#include "tmv/TMV_AuxProdXM.h"
  // Defines things like -m, x*m, m*x, x*(x*m), etc.
#undef PRODXM
#undef GENMATRIX

  //
  // Matrix + Scalar
  //

  template <class T, class T1, int N, StorageType S, IndexStyle I> 
  class SummX_1 : 
    public SmallMatrixComposite<T,N,N>
  {
    // m + x2
  public:
    inline SummX_1(
        T DEBUGPARAM(_x1), const SmallMatrix<T1,N,N,S,I>& _m, T _x2) :
      m(_m), x2(_x2) { TMVAssert(_x1 == T(1)); }
    inline const SmallMatrix<T1,N,N,S,I>& GetM() const 
    { return m; }
    inline T GetX2() const { return x2; }
    inline StorageType stor() const { return S; }
    inline void AssignTom(
        SmallMatrix<RealType(T),N,N,RowMajor,CStyle>& m0) const
    { TMVAssert(IsReal(T())); (m0=m) += REAL(x2); }
    inline void AssignTom(
        SmallMatrix<RealType(T),N,N,ColMajor,CStyle>& m0) const
    { TMVAssert(IsReal(T())); (m0=m) += REAL(x2); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),N,N,RowMajor,CStyle>& m0) const
    { (m0=m) += x2; }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),N,N,ColMajor,CStyle>& m0) const
    { (m0=m) += x2; }
    inline void AssignTom(
        SmallMatrix<RealType(T),N,N,RowMajor,FortranStyle>& m0) const
    { TMVAssert(IsReal(T())); (m0=m) += REAL(x2); }
    inline void AssignTom(
        SmallMatrix<RealType(T),N,N,ColMajor,FortranStyle>& m0) const
    { TMVAssert(IsReal(T())); (m0=m) += REAL(x2); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),N,N,RowMajor,FortranStyle>& m0) const
    { (m0=m) += x2; }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),N,N,ColMajor,FortranStyle>& m0) const
    { (m0=m) += x2; }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    {
      TMVAssert(IsReal(T()));
      if (m0.stepi() == 1 && m0.stepj() == N && !SameStorage(m0,m)) {
        DoCopy<N,N,S,ColMajor>(m.cptr(),m0.ptr());
        RealType(T)* m0p = m0.ptr();
        for(int i=0;i<N;i++) m0p[(N+1)*i] += REAL(x2);
      } else if (m0.stepj() == 1 && m0.stepi() == N && !SameStorage(m0,m)) {
        DoCopy<N,N,S,RowMajor>(m.cptr(),m0.ptr());
        RealType(T)* m0p = m0.ptr();
        for(int i=0;i<N;i++) m0p[(N+1)*i] += REAL(x2);
      } else {
        (m0=m) += REAL(x2); 
      }
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    {
      if (m0.stepi() == 1 && m0.stepj() == N && !SameStorage(m0,m)) {
        DoCopy<N,N,S,ColMajor>(m.cptr(),m0.ptr());
        ComplexType(T)* m0p = m0.ptr();
        for(int i=0;i<N;i++) m0p[(N+1)*i] += x2;
        if (m0.isconj()) m0.ConjugateSelf();
      } else if (m0.stepj() == 1 && m0.stepi() == N && !SameStorage(m0,m)) {
        DoCopy<N,N,S,RowMajor>(m.cptr(),m0.ptr());
        ComplexType(T)* m0p = m0.ptr();
        for(int i=0;i<N;i++) m0p[(N+1)*i] += x2;
        if (m0.isconj()) m0.ConjugateSelf();
      } else {
        (m0=m) += x2; 
      }
    }
  private:
    const SmallMatrix<T1,N,N,S,I>& m;
    const T x2;
  };

  template <class T, class T1, int N, StorageType S, IndexStyle I> 
  class SummX : 
    public SmallMatrixComposite<T,N,N>
  {
    // x1*m + x2
  public:
    inline SummX(T _x1, const SmallMatrix<T1,N,N,S,I>& _m, T _x2) :
      x1(_x1), m(_m), x2(_x2) {}
    inline T GetX1() const { return x1; }
    inline const SmallMatrix<T1,N,N,S,I>& GetM() const 
    { return m; }
    inline T GetX2() const { return x2; }
    inline StorageType stor() const { return S; }
    inline void AssignTom(
        SmallMatrix<RealType(T),N,N,RowMajor,CStyle>& m0) const
    { TMVAssert(IsReal(T())); (m0=x1*m) += REAL(x2); }
    inline void AssignTom(
        SmallMatrix<RealType(T),N,N,ColMajor,CStyle>& m0) const
    { TMVAssert(IsReal(T())); (m0=x1*m) += REAL(x2); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),N,N,RowMajor,CStyle>& m0) const
    { (m0=x1*m) += x2; }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),N,N,ColMajor,CStyle>& m0) const
    { (m0=x1*m) += x2; }
    inline void AssignTom(
        SmallMatrix<RealType(T),N,N,RowMajor,FortranStyle>& m0) const
    { TMVAssert(IsReal(T())); (m0=x1*m) += REAL(x2); }
    inline void AssignTom(
        SmallMatrix<RealType(T),N,N,ColMajor,FortranStyle>& m0) const
    { TMVAssert(IsReal(T())); (m0=x1*m) += REAL(x2); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),N,N,RowMajor,FortranStyle>& m0) const
    { (m0=x1*m) += x2; }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),N,N,ColMajor,FortranStyle>& m0) const
    { (m0=x1*m) += x2; }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    {
      TMVAssert(IsReal(T()));
      if (m0.stepi() == 1 && m0.stepj() == N && !SameStorage(m0,m)) {
        DoCopy<N,N,S,ColMajor>(m.cptr(),m0.ptr());
        RealType(T)* m0p = m0.ptr();
        MultXV<N*N>(x1,m0p);
        for(int i=0;i<N;i++) m0p[(N+1)*i] += REAL(x2);
      } else if (m0.stepj() == 1 && m0.stepi() == N && !SameStorage(m0,m)) {
        DoCopy<N,N,S,RowMajor>(m.cptr(),m0.ptr());
        RealType(T)* m0p = m0.ptr();
        MultXV<N*N>(x1,m0p);
        for(int i=0;i<N;i++) m0p[(N+1)*i] += REAL(x2);
      } else {
        (m0=x1*m) += REAL(x2); 
      }
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    {
      TMVAssert(IsReal(T()));
      if (m0.stepi() == 1 && m0.stepj() == N && !SameStorage(m0,m)) {
        DoCopy<N,N,S,ColMajor>(m.cptr(),m0.ptr());
        ComplexType(T)* m0p = m0.ptr();
        MultXV<N*N>(x1,m0p);
        for(int i=0;i<N;i++) m0p[(N+1)*i] += x2;
        if (m0.isconj()) m0.ConjugateSelf();
      } else if (m0.stepj() == 1 && m0.stepi() == N && !SameStorage(m0,m)) {
        DoCopy<N,N,S,RowMajor>(m.cptr(),m0.ptr());
        ComplexType(T)* m0p = m0.ptr();
        MultXV<N*N>(x1,m0p);
        for(int i=0;i<N;i++) m0p[(N+1)*i] += x2;
        if (m0.isconj()) m0.ConjugateSelf();
      } else {
        (m0=x1*m) += x2;
      }
    }
  private:
    const T x1;
    const SmallMatrix<T1,N,N,S,I>& m;
    const T x2;
  };

  // m+=x
  template <class T, int N, StorageType S, IndexStyle I>
  inline SmallMatrix<T,N,N,S,I>& operator+=(SmallMatrix<T,N,N,S,I>& m, T x) 
  { for(int i=0;i<N;i++) m.ref(i,i) += x; return m; }

  template <class T, int N, StorageType S, IndexStyle I>
  inline SmallMatrix<CT,N,N,S,I>& operator+=(SmallMatrix<CT,N,N,S,I>& m, T x) 
  { for(int i=0;i<N;i++) m.ref(i,i) += x; return m; }

  template <class T, int N, StorageType S, IndexStyle I>
  inline SmallMatrix<CT,N,N,S,I>& operator+=(SmallMatrix<CT,N,N,S,I>& m, CCT x) 
  { m += CT(x); return m; }

  template <class T, int N, StorageType S, IndexStyle I>
  inline SmallMatrix<CT,N,N,S,I>& operator+=(SmallMatrix<CT,N,N,S,I>& m, VCT x) 
  { m += CT(x); return m; }

  // m-=x
  template <class T, int N, StorageType S, IndexStyle I>
  inline SmallMatrix<T,N,N,S,I>& operator-=(SmallMatrix<T,N,N,S,I>& m, T x) 
  { m += (-x); return m; }

  template <class T, int N, StorageType S, IndexStyle I>
  inline SmallMatrix<CT,N,N,S,I>& operator-=(SmallMatrix<CT,N,N,S,I>& m, T x) 
  { m += (-x); return m; }

  template <class T, int N, StorageType S, IndexStyle I>
  inline SmallMatrix<CT,N,N,S,I>& operator-=(SmallMatrix<CT,N,N,S,I>& m, CCT x) 
  { m += CT(-x); return m; }

  template <class T, int N, StorageType S, IndexStyle I>
  inline SmallMatrix<CT,N,N,S,I>& operator-=(SmallMatrix<CT,N,N,S,I>& m, VCT x) 
  { m += CT(-x); return m; }

#define GENMATRIX SmallMatrix
#define PRODXM ProdXm
#define SUMMX SummX
#define SUMMX_1 SummX_1
#define X1 ,N,N,S,I
#define X2 ,N,S,I
#define Y ,int N, StorageType S, IndexStyle I
#include "tmv/TMV_AuxSumMX.h"
  // Defines things like m+x, x+m, x-m, m-x, x+x*m, x*(x+m), etc.
#undef SUMMX
#undef GENMATRIX
#undef PRODXM

  //
  // Vector ^ Vector (OuterProduct)
  //

  template <class T, class T1, class T2, int M, int N, IndexStyle I1, IndexStyle I2>
  class OProdvv_1 : 
    public SmallMatrixComposite<T,M,N>
  {
  public:
    inline OProdvv_1(
        const T DEBUGPARAM(_x), const SmallVector<T1,M,I1>& _v1,
        const SmallVector<T2,N,I2>& _v2) : 
      v1(_v1), v2(_v2) { TMVAssert(_x == T(1)); }
    inline StorageType stor() const { return ColMajor; }
    inline const SmallVector<T1,M,I1>& GetV1() const { return v1; }
    inline const SmallVector<T2,N,I2>& GetV2() const { return v2; }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,CStyle>& m0) const
    { 
      TMVAssert(IsReal(T()));
      Rank1Update_1<M,N,RowMajor>(v1.cptr(),v2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,CStyle>& m0) const
    { 
      TMVAssert(IsReal(T()));
      Rank1Update_1<M,N,ColMajor>(v1.cptr(),v2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,CStyle>& m0) const
    { Rank1Update_1<M,N,RowMajor>(v1.cptr(),v2.cptr(),m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,CStyle>& m0) const
    { Rank1Update_1<M,N,ColMajor>(v1.cptr(),v2.cptr(),m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,FortranStyle>& m0) const
    { 
      TMVAssert(IsReal(T()));
      Rank1Update_1<M,N,RowMajor>(v1.cptr(),v2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,FortranStyle>& m0) const
    { 
      TMVAssert(IsReal(T()));
      Rank1Update_1<M,N,ColMajor>(v1.cptr(),v2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,FortranStyle>& m0) const
    { Rank1Update_1<M,N,RowMajor>(v1.cptr(),v2.cptr(),m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,FortranStyle>& m0) const
    { Rank1Update_1<M,N,ColMajor>(v1.cptr(),v2.cptr(),m0.ptr()); }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    {
      TMVAssert(IsReal(T()));
      if (m0.stepi() == 1 && m0.stepj() == M) {
        Rank1Update_1<M,N,ColMajor>(v1.cptr(),v2.cptr(),m0.ptr());
      } else if (m0.stepj() == 1 && m0.stepi() == N) {
        Rank1Update_1<M,N,ColMajor>(v1.cptr(),v2.cptr(),m0.ptr());
      } else {
        Rank1Update<false>(T(1),v1.View(),v2.View(),m0.View()); 
      }
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    {
      if (m0.stepi() == 1 && m0.stepj() == M && !m0.isconj()) {
        Rank1Update_1<M,N,ColMajor>(v1.cptr(),v2.cptr(),m0.ptr());
      } else if (m0.stepj() == 1 && m0.stepi() == N && !m0.isconj()) {
        Rank1Update_1<M,N,ColMajor>(v1.cptr(),v2.cptr(),m0.ptr());
      } else {
        Rank1Update<false>(T(1),v1.View(),v2.View(),m0.View()); 
      }
    }
  private:
    const SmallVector<T1,M,I1>& v1;
    const SmallVector<T2,N,I2>& v2;
  };

  template <class T, class T1, class T2, int M, int N, IndexStyle I1, IndexStyle I2>
  class OProdvv : 
    public SmallMatrixComposite<T,M,N>
  {
  public:
    inline OProdvv(const T _x, const SmallVector<T1,M,I1>& _v1,
        const SmallVector<T2,N,I2>& _v2) :
      x(_x), v1(_v1), v2(_v2) {}
    inline StorageType stor() const { return ColMajor; }
    inline T GetX() const { return x; }
    inline const SmallVector<T1,M,I1>& GetV1() const { return v1; }
    inline const SmallVector<T2,N,I2>& GetV2() const { return v2; }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,CStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      Rank1Update<M,N,RowMajor>(x,v1.cptr(),v2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,CStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      Rank1Update<M,N,ColMajor>(x,v1.cptr(),v2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,CStyle>& m0) const
    { Rank1Update<M,N,RowMajor>(x,v1.cptr(),v2.cptr(),m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,CStyle>& m0) const
    { Rank1Update<M,N,ColMajor>(x,v1.cptr(),v2.cptr(),m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,FortranStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      Rank1Update<M,N,RowMajor>(x,v1.cptr(),v2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,FortranStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      Rank1Update<M,N,ColMajor>(x,v1.cptr(),v2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,FortranStyle>& m0) const
    { Rank1Update<M,N,RowMajor>(x,v1.cptr(),v2.cptr(),m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,FortranStyle>& m0) const
    { Rank1Update<M,N,ColMajor>(x,v1.cptr(),v2.cptr(),m0.ptr()); }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    {
      TMVAssert(IsReal(T()));
      if (m0.stepi() == 1 && m0.stepj() == M) {
        Rank1Update<M,N,ColMajor>(x,v1.cptr(),v2.cptr(),m0.ptr());
      } else if (m0.stepj() == 1 && m0.stepi() == N) {
        Rank1Update<M,N,ColMajor>(x,v1.cptr(),v2.cptr(),m0.ptr());
      } else {
        Rank1Update<false>(x,v1.View(),v2.View(),m0.View()); 
      }
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    { 
      if (m0.stepi() == 1 && m0.stepj() == M && !m0.isconj()) {
        Rank1Update<M,N,ColMajor>(x,v1.cptr(),v2.cptr(),m0.ptr());
      } else if (m0.stepj() == 1 && m0.stepi() == N && !m0.isconj()) {
        Rank1Update<M,N,ColMajor>(x,v1.cptr(),v2.cptr(),m0.ptr());
      } else {
        Rank1Update<false>(x,v1.View(),v2.View(),m0.View()); 
      }
    }
  private:
    T x;
    const SmallVector<T1,M,I1>& v1;
    const SmallVector<T2,N,I2>& v2;
  };

  template <class T, class T1, class T2, int N, IndexStyle I>
  class OProdVv : 
    public MatrixComposite<T>
  {
  public:
    inline OProdVv(const T _x, const GenVector<T1>& _v1,
        const SmallVector<T2,N,I>& _v2) :
      x(_x), v1(_v1), v2(_v2) {}
    inline size_t colsize() const { return v1.size(); }
    inline size_t rowsize() const { return N; }
    inline StorageType stor() const { return ColMajor; }
    inline T GetX() const { return x; }
    inline const GenVector<T1>& GetV1() const { return v1; }
    inline const SmallVector<T2,N,I>& GetV2() const { return v2; }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    {
      TMVAssert(IsReal(T()));
      Rank1Update<false>(x,v1,v2.View(),m0.View()); 
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    { Rank1Update<false>(x,v1,v2.View(),m0.View()); }
  private:
    T x;
    const GenVector<T1>& v1;
    const SmallVector<T2,N,I>& v2;
  };

  template <class T, class T1, class T2, int M, IndexStyle I>
  class OProdvV : 
    public MatrixComposite<T>
  {
  public:
    inline OProdvV(const T _x, const SmallVector<T1,M,I>& _v1,
        const GenVector<T2>& _v2) :
      x(_x), v1(_v1), v2(_v2) {}
    inline size_t colsize() const { return M; }
    inline size_t rowsize() const { return v2.size(); }
    inline StorageType stor() const { return ColMajor; }
    inline T GetX() const { return x; }
    inline const SmallVector<T1,M,I>& GetV1() const { return v1; }
    inline const GenVector<T2>& GetV2() const { return v2; }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    {
      TMVAssert(IsReal(T()));
      Rank1Update<false>(x,v1.View(),v2,m0.View()); 
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    { Rank1Update<false>(x,v1.View(),v2,m0.View()); }
  private:
    T x;
    const SmallVector<T1,M,I>& v1;
    const GenVector<T2>& v2;
  };

  // m+=(v^v)
  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I, IndexStyle I1, IndexStyle I2>
  inline SmallMatrix<T,M,N,S,I>& operator+=(
      SmallMatrix<T,M,N,S,I>& m0, const OProdvv_1<T,T1,T2,M,N,I1,I2>& opvv)
  { 
    AddRank1Update_1<M,N,S>(opvv.GetV1().cptr(), opvv.GetV2().cptr(), 
        m0.ptr());
    return m0;
  }

  template <class T, int M, int N, StorageType S, IndexStyle I, IndexStyle I1, IndexStyle I2>
  inline SmallMatrix<CT,M,N,S,I>& operator+=(
      SmallMatrix<CT,M,N,S,I>& m0, const OProdvv_1<T,T,T,M,N,I1,I2>& opvv)
  { 
    AddRank1Update_1<M,N,S>(opvv.GetV1().cptr(), opvv.GetV2().cptr(), 
        m0.ptr());
    return m0;
  }

  template <class T, class T1, class T2, int M, int N, IndexStyle I1, IndexStyle I2>
  inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m0, const OProdvv_1<T,T1,T2,M,N,I1,I2>& opvv)
  {
    TMVAssert(m0.colsize() == opvv.colsize());
    TMVAssert(m0.rowsize() == opvv.rowsize());
    Rank1Update<true>(T(1),opvv.GetV1().View(), opvv.GetV2().View(), m0);
    return m0; 
  }

  template <class T, int M, int N, IndexStyle I1, IndexStyle I2> 
  inline const MatrixView<CT>& operator+=(
      const MatrixView<CT>& m0, const OProdvv_1<T,T,T,M,N,I1,I2>& opvv)
  { 
    TMVAssert(m0.colsize() == opvv.colsize());
    TMVAssert(m0.rowsize() == opvv.rowsize());
    Rank1Update<true>(CT(1),opvv.GetV1().View(), opvv.GetV2().View(), m0); 
    return m0; 
  }

  // m-=(v^v)
  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I, IndexStyle I1, IndexStyle I2>
  inline SmallMatrix<T,M,N,S,I>& operator-=(
      SmallMatrix<T,M,N,S,I>& m0, const OProdvv_1<T,T1,T2,M,N,I1,I2>& opvv)
  { 
    AddRank1Update_m1<M,N,S>(opvv.GetV1().cptr(), opvv.GetV2().cptr(), 
        m0.ptr());
    return m0;
  }

  template <class T, int M, int N, StorageType S, IndexStyle I, IndexStyle I1, IndexStyle I2>
  inline SmallMatrix<CT,M,N,S,I>& operator-=(
      SmallMatrix<CT,M,N,S,I>& m0, const OProdvv_1<T,T,T,M,N,I1,I2>& opvv)
  { 
    AddRank1Update_m1<M,N,S>(opvv.GetV1().cptr(), opvv.GetV2().cptr(), 
        m0.ptr());
    return m0;
  }

  template <class T, class T1, class T2, int M, int N, IndexStyle I1, IndexStyle I2>
  inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m0, const OProdvv_1<T,T1,T2,M,N,I1,I2>& opvv)
  { 
    TMVAssert(m0.colsize() == opvv.colsize());
    TMVAssert(m0.rowsize() == opvv.rowsize());
    Rank1Update<true>(T(-1),opvv.GetV1().View(), opvv.GetV2().View(), m0);
    return m0; 
  }

  template <class T, int M, int N, IndexStyle I1, IndexStyle I2> 
  inline const MatrixView<CT>& operator-=(
      const MatrixView<CT>& m0, const OProdvv_1<T,T,T,M,N,I1,I2>& opvv)
  {
    TMVAssert(m0.colsize() == opvv.colsize());
    TMVAssert(m0.rowsize() == opvv.rowsize());
    Rank1Update<true>(CT(-1), opvv.GetV1().View(), opvv.GetV2().View(), m0);
    return m0; 
  }

  // m+=(x*v^v)
  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I, IndexStyle I1, IndexStyle I2>
  inline SmallMatrix<T,M,N,S,I>& operator+=(
      SmallMatrix<T,M,N,S,I>& m0, const OProdvv<T,T1,T2,M,N,I1,I2>& opvv)
  {
    AddRank1Update<M,N,S>(opvv.GetX(), opvv.GetV1().cptr(), 
        opvv.GetV2().cptr(), m0.ptr());
    return m0;
  }

  template <class T, int M, int N, StorageType S, IndexStyle I, IndexStyle I1, IndexStyle I2>
  inline SmallMatrix<CT,M,N,S,I>& operator+=(
      SmallMatrix<CT,M,N,S,I>& m0, const OProdvv<T,T,T,M,N,I1,I2>& opvv)
  {
    AddRank1Update<M,N,S>(opvv.GetX(), opvv.GetV1().cptr(),
        opvv.GetV2().cptr(), m0.ptr());
    return m0;
  }

  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I>
  inline SmallMatrix<T,M,N,S,I>& operator+=(
      SmallMatrix<T,M,N,S,I>& m0, const OProdVV<T,T1,T2>& opvv)
  { 
    TMVAssert(m0.colsize() == opvv.colsize());
    TMVAssert(m0.rowsize() == opvv.rowsize());
    Rank1Update<true>(opvv.GetX(), opvv.GetV1(), opvv.GetV2(), m0.View());
    return m0;
  }

  template <class T, int M, int N, StorageType S, IndexStyle I>
  inline SmallMatrix<CT,M,N,S,I>& operator+=(
      SmallMatrix<CT,M,N,S,I>& m0, const OProdVV<T,T,T>& opvv)
  { 
    TMVAssert(m0.colsize() == opvv.colsize());
    TMVAssert(m0.rowsize() == opvv.rowsize());
    Rank1Update<true>(opvv.GetX(), opvv.GetV1(), opvv.GetV2(), m0.View());
    return m0;
  }

  template <class T, class T1, class T2, int M, int N, IndexStyle I1, IndexStyle I2>
  inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m0, const OProdvv<T,T1,T2,M,N,I1,I2>& opvv)
  { 
    TMVAssert(m0.colsize() == opvv.colsize());
    TMVAssert(m0.rowsize() == opvv.rowsize());
    Rank1Update<true>(opvv.GetX(), opvv.GetV1().View(), 
        opvv.GetV2().View(), m0);
    return m0; 
  }

  template <class T, int M, int N, IndexStyle I1, IndexStyle I2> 
  inline const MatrixView<CT>& operator+=(
      const MatrixView<CT>& m0, const OProdvv<T,T,T,M,N,I1,I2>& opvv)
  {
    TMVAssert(m0.colsize() == opvv.colsize());
    TMVAssert(m0.rowsize() == opvv.rowsize());
    Rank1Update<true>(opvv.GetX(), opvv.GetV1().View(),
        opvv.GetV2().View(), m0);
    return m0; 
  }

  // m-=(x*v^v)
  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I, IndexStyle I1, IndexStyle I2>
  inline SmallMatrix<T,M,N,S,I>& operator-=(
      SmallMatrix<T,M,N,S,I>& m0, const OProdvv<T,T1,T2,M,N,I1,I2>& opvv)
  {
    AddRank1Update<M,N,S>(-opvv.GetX(), opvv.GetV1().cptr(), 
        opvv.GetV2().cptr(), m0.ptr());
    return m0;
  }

  template <class T, int M, int N, StorageType S, IndexStyle I, IndexStyle I1, IndexStyle I2>
  inline SmallMatrix<CT,M,N,S,I>& operator-=(
      SmallMatrix<CT,M,N,S,I>& m0, const OProdvv<T,T,T,M,N,I1,I2>& opvv)
  { 
    AddRank1Update<M,N,S>(-opvv.GetX(), opvv.GetV1().cptr(), 
        opvv.GetV2().cptr(), m0.ptr());
    return m0;
  }

  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I>
  inline SmallMatrix<T,M,N,S,I>& operator-=(
      SmallMatrix<T,M,N,S,I>& m0, const OProdVV<T,T1,T2>& opvv)
  { 
    TMVAssert(m0.colsize() == opvv.colsize());
    TMVAssert(m0.rowsize() == opvv.rowsize());
    Rank1Update<true>(-opvv.GetX(), opvv.GetV1(), opvv.GetV2(), m0.View());
    return m0;
  }

  template <class T, int M, int N, StorageType S, IndexStyle I>
  inline SmallMatrix<CT,M,N,S,I>& operator-=(
      SmallMatrix<CT,M,N,S,I>& m0, const OProdVV<T,T,T>& opvv)
  { 
    TMVAssert(m0.colsize() == opvv.colsize());
    TMVAssert(m0.rowsize() == opvv.rowsize());
    Rank1Update<true>(-opvv.GetX(), opvv.GetV1(), opvv.GetV2(), m0.View());
    return m0;
  }

  template <class T, class T1, class T2, int M, int N, IndexStyle I1, IndexStyle I2>
  inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m0, const OProdvv<T,T1,T2,M,N,I1,I2>& opvv)
  { 
    TMVAssert(m0.colsize() == opvv.colsize());
    TMVAssert(m0.rowsize() == opvv.rowsize());
    Rank1Update<true>(-opvv.GetX(), opvv.GetV1().View(),
        opvv.GetV2().View(), m0);
    return m0; 
  }

  template <class T, int M, int N, IndexStyle I1, IndexStyle I2> 
  inline const MatrixView<CT>& operator-=(
      const MatrixView<CT>& m0, const OProdvv<T,T,T,M,N,I1,I2>& opvv)
  {
    TMVAssert(m0.colsize() == opvv.colsize());
    TMVAssert(m0.rowsize() == opvv.rowsize());
    Rank1Update<true>(-opvv.GetX(), opvv.GetV1().View(),
        opvv.GetV2().View(), m0);
    return m0; 
  }

#define PRODMM OProdvv
#define PRODMM_1 OProdvv_1
#define GENMATRIX1 SmallVector
#define GENMATRIX2 SmallVector
#define PRODXM1 ProdXv
#define PRODXM2 ProdXv
#define X1 ,M,I1
#define X2 ,N,I2
#define X3 ,M,N,I1,I2
#define Y ,int M,int N, IndexStyle I1, IndexStyle I2
#define OP operator^
#define GETM1 .GetV()
#define GETM2 .GetV()
#include "tmv/TMV_AuxProdMM.h"
#define PRODMM_1 OProdvv_1
#define X3 ,M,N,I1,I2
#define Y ,int M,int N, IndexStyle I1, IndexStyle I2
#define GETM1 .GetV1()
#define GETM2 .GetV2()
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

#define PRODMM OProdVv
#define GENMATRIX1 GenVector
#define GENMATRIX2 SmallVector
#define PRODXM1 ProdXV
#define PRODXM2 ProdXv
#define X2 ,N,I
#define X3 ,N,I
#define Y ,int N, IndexStyle I
#define OP operator^
#define GETM1 .GetV()
#define GETM2 .GetV()
#include "tmv/TMV_AuxProdMM.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

#define PRODMM OProdvV
#define GENMATRIX1 SmallVector
#define GENMATRIX2 GenVector
#define PRODXM1 ProdXv
#define PRODXM2 ProdXV
#define X1 ,M,I
#define X3 ,M,I
#define Y ,int M, IndexStyle I
#define OP operator^
#define GETM1 .GetV()
#define GETM2 .GetV()
#include "tmv/TMV_AuxProdMM.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2


  //
  // Matrix + Matrix
  //

  template <class T, class T1, class T2, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
  class Summm_1_1 : 
    public SmallMatrixComposite<T,M,N> 
  {
  public:
    inline Summm_1_1(
        const T DEBUGPARAM(_x1), const SmallMatrix<T1,M,N,S1,I1>& _m1, 
        const T DEBUGPARAM(_x2), const SmallMatrix<T2,M,N,S2,I2>& _m2) :
      m1(_m1), m2(_m2) { TMVAssert(_x1 == T(1) && _x2 == T(1)); }
    inline StorageType stor() const { return S1; }
    inline const SmallMatrix<T1,M,N,S1,I1>& GetM1() const 
    { return m1; }
    inline const SmallMatrix<T2,M,N,S2,I2>& GetM2() const 
    { return m2; }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,CStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0 = m1;
      AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,CStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0 = m1;
      AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,CStyle>& m0) const
    { 
      m0 = m1;
      AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,CStyle>& m0) const
    { 
      m0 = m1;
      AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,FortranStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0 = m1;
      AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,FortranStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0 = m1;
      AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,FortranStyle>& m0) const
    { 
      m0 = m1;
      AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,FortranStyle>& m0) const
    { 
      m0 = m1;
      AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    {
      TMVAssert(IsReal(T()));
      if (m0.stepi() == 1 && m0.stepj() == M && 
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        DoCopy<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
        AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
      } else if (m0.stepj() == 1 && m0.stepi() == N &&
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        DoCopy<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
        AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
      } else {
        AddMM(T(1),m1.View(),T(1),m2.View(),m0); 
      }
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    {
      if (m0.stepi() == 1 && m0.stepj() == M && 
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        DoCopy<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
        AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
        if (m0.isconj()) m0.ConjugateSelf();
      } else if (m0.stepj() == 1 && m0.stepi() == N &&
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        DoCopy<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
        AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
        if (m0.isconj()) m0.ConjugateSelf();
      } else {
        AddMM(T(1),m1.View(),T(1),m2.View(),m0); 
      }
    }
  private:
    const SmallMatrix<T1,M,N,S1,I1>& m1;
    const SmallMatrix<T2,M,N,S2,I2>& m2;
  };

  template <class T, class T1, class T2, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
  class Summm_1_m1 : 
    public SmallMatrixComposite<T,M,N> 
  {
  public:
    inline Summm_1_m1(
        const T DEBUGPARAM(_x1), const SmallMatrix<T1,M,N,S1,I1>& _m1, 
        const T DEBUGPARAM(_x2), const SmallMatrix<T2,M,N,S2,I2>& _m2) :
      m1(_m1), m2(_m2) { TMVAssert(_x1 == T(1) && _x2 == T(-1)); }
    inline StorageType stor() const { return S1; }
    inline const SmallMatrix<T1,M,N,S1,I1>& GetM1() const { return m1; }
    inline const SmallMatrix<T2,M,N,S2,I2>& GetM2() const { return m2; }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,CStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0 = m1;
      AddMM_m1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,CStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0 = m1;
      AddMM_m1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,CStyle>& m0) const
    { 
      m0 = m1;
      AddMM_m1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,CStyle>& m0) const
    { 
      m0 = m1;
      AddMM_m1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,FortranStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0 = m1;
      AddMM_m1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,FortranStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0 = m1;
      AddMM_m1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,FortranStyle>& m0) const
    { 
      m0 = m1;
      AddMM_m1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,FortranStyle>& m0) const
    { 
      m0 = m1;
      AddMM_m1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    {
      TMVAssert(IsReal(T()));
      if (m0.stepi() == 1 && m0.stepj() == M && 
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        DoCopy<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
        AddMM_m1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
      } else if (m0.stepj() == 1 && m0.stepi() == N &&
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        DoCopy<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
        AddMM_m1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
      } else {
        AddMM(T(1),m1.View(),T(-1),m2.View(),m0); 
      }
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    {
      if (m0.stepi() == 1 && m0.stepj() == M && 
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        DoCopy<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
        AddMM_m1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
        if (m0.isconj()) m0.ConjugateSelf();
      } else if (m0.stepj() == 1 && m0.stepi() == N &&
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        DoCopy<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
        AddMM_m1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
        if (m0.isconj()) m0.ConjugateSelf();
      } else {
        AddMM(T(1),m1.View(),T(-1),m2.View(),m0); 
      }
    }
  private:
    const SmallMatrix<T1,M,N,S1,I1>& m1;
    const SmallMatrix<T2,M,N,S2,I2>& m2;
  };

  template <class T, class T1, class T2, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
  class Summm_1_x : 
    public SmallMatrixComposite<T,M,N> 
  {
  public:
    inline Summm_1_x(
        const T DEBUGPARAM(_x1), const SmallMatrix<T1,M,N,S1,I1>& _m1, 
        const T _x2, const SmallMatrix<T2,M,N,S2,I2>& _m2) :
      m1(_m1), x2(_x2), m2(_m2) { TMVAssert(_x1 == T(1)); }
    inline StorageType stor() const { return S1; }
    inline const SmallMatrix<T1,M,N,S1,I1>& GetM1() const { return m1; }
    inline T GetX2() const { return x2; }
    inline const SmallMatrix<T2,M,N,S2,I2>& GetM2() const { return m2; }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,CStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0 = x2*m2;
      AddMM_1<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,CStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0 = x2*m2;
      AddMM_1<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,CStyle>& m0) const
    { 
      m0 = x2*m2;
      AddMM_1<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,CStyle>& m0) const
    { 
      m0 = x2*m2;
      AddMM_1<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,FortranStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0 = x2*m2;
      AddMM_1<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,FortranStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0 = x2*m2;
      AddMM_1<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,FortranStyle>& m0) const
    { 
      m0 = x2*m2;
      AddMM_1<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,FortranStyle>& m0) const
    { 
      m0 = x2*m2;
      AddMM_1<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
    }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    {
      TMVAssert(IsReal(T()));
      if (m0.stepi() == 1 && m0.stepj() == M && 
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        DoCopy<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
        MultXV<M*N>(x2,m0.ptr());
        AddMM_1<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
      } else if (m0.stepj() == 1 && m0.stepi() == N &&
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        DoCopy<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
        MultXV<M*N>(x2,m0.ptr());
        AddMM_1<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
      } else {
        AddMM(T(1),m1.View(),x2,m2.View(),m0); 
      }
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    {
      if (m0.stepi() == 1 && m0.stepj() == M && 
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        DoCopy<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
        MultXV<M*N>(x2,m0.ptr());
        AddMM_1<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
        if (m0.isconj()) m0.ConjugateSelf();
      } else if (m0.stepj() == 1 && m0.stepi() == N &&
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        DoCopy<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
        MultXV<M*N>(x2,m0.ptr());
        AddMM_1<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
        if (m0.isconj()) m0.ConjugateSelf();
      } else {
        AddMM(T(1),m1.View(),x2,m2.View(),m0); 
      }
    }
  private:
    const SmallMatrix<T1,M,N,S1,I1>& m1;
    const T x2;
    const SmallMatrix<T2,M,N,S2,I2>& m2;
  };

  template <class T, class T1, class T2, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
  class Summm_x_1 : 
    public SmallMatrixComposite<T,M,N> 
  {
  public:
    inline Summm_x_1(
        const T _x1, const SmallMatrix<T1,M,N,S1,I1>& _m1, 
        const T DEBUGPARAM(_x2), const SmallMatrix<T2,M,N,S2,I2>& _m2) :
      x1(_x1), m1(_m1), m2(_m2) { TMVAssert(_x2 == T(1)); }
    inline StorageType stor() const { return S1; }
    inline T GetX1() const { return x1; }
    inline const SmallMatrix<T1,M,N,S1,I1>& GetM1() const { return m1; }
    inline const SmallMatrix<T2,M,N,S2,I2>& GetM2() const { return m2; }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,CStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0 = x1*m1;
      AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,CStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0 = x1*m1;
      AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,CStyle>& m0) const
    { 
      m0 = x1*m1;
      AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,CStyle>& m0) const
    { 
      m0 = x1*m1;
      AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,FortranStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0 = x1*m1;
      AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,FortranStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0 = x1*m1;
      AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,FortranStyle>& m0) const
    { 
      m0 = x1*m1;
      AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,FortranStyle>& m0) const
    { 
      m0 = x1*m1;
      AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    {
      TMVAssert(IsReal(T()));
      if (m0.stepi() == 1 && m0.stepj() == M && 
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        DoCopy<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
        MultXV<M*N>(x1,m0.ptr());
        AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
      } else if (m0.stepj() == 1 && m0.stepi() == N &&
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        DoCopy<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
        MultXV<M*N>(x1,m0.ptr());
        AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
      } else {
        AddMM(x1,m1.View(),T(1),m2.View(),m0); 
      }
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    {
      if (m0.stepi() == 1 && m0.stepj() == M && 
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        DoCopy<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
        MultXV<M*N>(x1,m0.ptr());
        AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
        if (m0.isconj()) m0.ConjugateSelf();
      } else if (m0.stepj() == 1 && m0.stepi() == N &&
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        DoCopy<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
        MultXV<M*N>(x1,m0.ptr());
        AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
        if (m0.isconj()) m0.ConjugateSelf();
      } else {
        AddMM(x1,m1.View(),T(1),m2.View(),m0); 
      }
    }
  private:
    const T x1;
    const SmallMatrix<T1,M,N,S1,I1>& m1;
    const SmallMatrix<T2,M,N,S2,I2>& m2;
  };

  template <class T, class T1, class T2, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
  class Summm_x_m1 : 
    public SmallMatrixComposite<T,M,N> 
  {
  public:
    inline Summm_x_m1(
        const T _x1, const SmallMatrix<T1,M,N,S1,I1>& _m1, 
        const T DEBUGPARAM(_x2), const SmallMatrix<T2,M,N,S2,I2>& _m2) :
      x1(_x1), m1(_m1), m2(_m2) { TMVAssert(_x2 == T(-1)); }
    inline StorageType stor() const { return S1; }
    inline T GetX1() const { return x1; }
    inline const SmallMatrix<T1,M,N,S1,I1>& GetM1() const { return m1; }
    inline const SmallMatrix<T2,M,N,S2,I2>& GetM2() const { return m2; }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,CStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0 = x1*m1;
      AddMM_m1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,CStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0 = x1*m1;
      AddMM_m1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,CStyle>& m0) const
    { 
      m0 = x1*m1;
      AddMM_m1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,CStyle>& m0) const
    { 
      m0 = x1*m1;
      AddMM_m1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,FortranStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0 = x1*m1;
      AddMM_m1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,FortranStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0 = x1*m1;
      AddMM_m1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,FortranStyle>& m0) const
    { 
      m0 = x1*m1;
      AddMM_m1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,FortranStyle>& m0) const
    { 
      m0 = x1*m1;
      AddMM_m1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
    }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    {
      TMVAssert(IsReal(T()));
      if (m0.stepi() == 1 && m0.stepj() == M && 
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        DoCopy<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
        MultXV<M*N>(x1,m0.ptr());
        AddMM_m1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
      } else if (m0.stepj() == 1 && m0.stepi() == N &&
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        DoCopy<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
        MultXV<M*N>(x1,m0.ptr());
        AddMM_m1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
      } else {
        AddMM(x1,m1.View(),T(1),m2.View(),m0); 
      }
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    {
      if (m0.stepi() == 1 && m0.stepj() == M && 
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        DoCopy<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
        MultXV<M*N>(x1,m0.ptr());
        AddMM_m1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
        if (m0.isconj()) m0.ConjugateSelf();
      } else if (m0.stepj() == 1 && m0.stepi() == N &&
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        DoCopy<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
        MultXV<M*N>(x1,m0.ptr());
        AddMM_m1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
        if (m0.isconj()) m0.ConjugateSelf();
      } else {
        AddMM(x1,m1.View(),T(1),m2.View(),m0); 
      }
    }
  private:
    const T x1;
    const SmallMatrix<T1,M,N,S1,I1>& m1;
    const SmallMatrix<T2,M,N,S2,I2>& m2;
  };

  template <class T, class T1, class T2, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
  class Summm : 
    public SmallMatrixComposite<T,M,N> 
  {
  public:
    inline Summm(const T _x1, const SmallMatrix<T1,M,N,S1,I1>& _m1, 
        const T _x2, const SmallMatrix<T2,M,N,S2,I2>& _m2) :
      x1(_x1), m1(_m1), x2(_x2), m2(_m2) {}
    inline StorageType stor() const { return S1; }
    inline T GetX1() const { return x1; }
    inline const SmallMatrix<T1,M,N,S1,I1>& GetM1() const { return m1; }
    inline T GetX2() const { return x2; }
    inline const SmallMatrix<T2,M,N,S2,I2>& GetM2() const { return m2; }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,CStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0 = x1*m1;
      AddMM<M,N,S2,RowMajor>(x2,m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,CStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0 = x1*m1;
      AddMM<M,N,S2,ColMajor>(x2,m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,CStyle>& m0) const
    { 
      m0 = x1*m1;
      AddMM<M,N,S2,RowMajor>(x2,m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,CStyle>& m0) const
    { 
      m0 = x1*m1;
      AddMM<M,N,S2,ColMajor>(x2,m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,FortranStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0 = x1*m1;
      AddMM<M,N,S2,RowMajor>(x2,m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,FortranStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      m0 = x1*m1;
      AddMM<M,N,S2,ColMajor>(x2,m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,FortranStyle>& m0) const
    { 
      m0 = x1*m1;
      AddMM<M,N,S2,RowMajor>(x2,m2.cptr(),m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,FortranStyle>& m0) const
    { 
      m0 = x1*m1;
      AddMM<M,N,S2,ColMajor>(x2,m2.cptr(),m0.ptr());
    }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    {
      TMVAssert(IsReal(T()));
      if (m0.stepi() == 1 && m0.stepj() == M && 
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        DoCopy<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
        MultXV<M*N>(x1,m0.ptr());
        AddMM<M,N,S2,ColMajor>(x2,m2.cptr(),m0.ptr());
      } else if (m0.stepj() == 1 && m0.stepi() == N &&
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        DoCopy<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
        MultXV<M*N>(x1,m0.ptr());
        AddMM<M,N,S2,RowMajor>(x2,m2.cptr(),m0.ptr());
      } else {
        AddMM(x1,m1.View(),x2,m2.View(),m0); 
      }
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    {
      if (m0.stepi() == 1 && m0.stepj() == M && 
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        DoCopy<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
        MultXV<M*N>(x1,m0.ptr());
        AddMM<M,N,S2,ColMajor>(x2,m2.cptr(),m0.ptr());
        if (m0.isconj()) m0.ConjugateSelf();
      } else if (m0.stepj() == 1 && m0.stepi() == N &&
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        DoCopy<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
        MultXV<M*N>(x1,m0.ptr());
        AddMM<M,N,S2,RowMajor>(x2,m2.cptr(),m0.ptr());
        if (m0.isconj()) m0.ConjugateSelf();
      } else {
        AddMM(x1,m1.View(),x2,m2.View(),m0); 
      }
    }
  private:
    const T x1;
    const SmallMatrix<T1,M,N,S1,I1>& m1;
    const T x2;
    const SmallMatrix<T2,M,N,S2,I2>& m2;
  };

  template <class T, class T1, class T2, int M, int N, StorageType S1, IndexStyle I1> 
  class SummM : 
    public MatrixComposite<T> 
  {
  public:
    inline SummM(const T _x1, const SmallMatrix<T1,M,N,S1,I1>& _m1, 
        const T _x2, const GenMatrix<T2>& _m2) :
      x1(_x1), m1(_m1), x2(_x2), m2(_m2) {}
    inline size_t colsize() const { return M; }
    inline size_t rowsize() const { return N; }
    inline StorageType stor() const { return S1; }
    inline T GetX1() const { return x1; }
    inline const SmallMatrix<T1,M,N,S1,I1>& GetM1() const { return m1; }
    inline T GetX2() const { return x2; }
    inline const GenMatrix<T2>& GetM2() const { return m2; }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    {
      TMVAssert(IsReal(T()));
      if (m2.stepi() == 1 && m2.stepj() == M && 
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        if (m0.stepi() == 1 && m0.stepj() == M) {
          DoCopy<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
          if (x1 != T(1)) MultXV<M*N>(x1,m0.ptr());
          if (x2 == T(1))
            AddMM_1<M,N,ColMajor,ColMajor>(m2.cptr(),m0.ptr());
          else 
            AddMM<M,N,ColMajor,ColMajor>(x2,m2.cptr(),m0.ptr());
        } else if (m0.stepj() == 1 && m0.stepi() == N) {
          DoCopy<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
          if (x1 != T(1)) MultXV<M*N>(x1,m0.ptr());
          if (x2 == T(1))
            AddMM_1<M,N,ColMajor,RowMajor>(m2.cptr(),m0.ptr());
          else 
            AddMM<M,N,ColMajor,RowMajor>(x2,m2.cptr(),m0.ptr());
        } else {
          AddMM(x1,m1.View(),x2,m2.View(),m0); 
        }
      } else if (m2.stepj() == 1 && m2.stepi() == N && 
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        if (m0.stepi() == 1 && m0.stepj() == M) {
          DoCopy<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
          if (x1 != T(1)) MultXV<M*N>(x1,m0.ptr());
          if (x2 == T(1))
            AddMM_1<M,N,RowMajor,ColMajor>(m2.cptr(),m0.ptr());
          else 
            AddMM<M,N,RowMajor,ColMajor>(x2,m2.cptr(),m0.ptr());
        } else if (m0.stepj() == 1 && m0.stepi() == N) {
          DoCopy<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
          if (x1 != T(1)) MultXV<M*N>(x1,m0.ptr());
          if (x2 == T(1))
            AddMM_1<M,N,RowMajor,RowMajor>(m2.cptr(),m0.ptr());
          else 
            AddMM<M,N,RowMajor,RowMajor>(x2,m2.cptr(),m0.ptr());
        } else {
          AddMM(x1,m1.View(),x2,m2.View(),m0); 
        }
      } else {
        AddMM(x1,m1.View(),x2,m2.View(),m0); 
      }
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    {
      if (m2.stepi() == 1 && m2.stepj() == M &&
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        if (m0.stepi() == 1 && m0.stepj() == M) {
          DoCopy<M,N,ColMajor,ColMajor>(m2.cptr(),m0.ptr());
          if (m2.isconj()) m0.ConjugateSelf();
          if (x2 != T(1)) MultXV<M*N>(x2,m0.ptr());
          if (x1 == T(1))
            AddMM_1<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
          else 
            AddMM<M,N,S1,ColMajor>(x1,m1.cptr(),m0.ptr());
          if (m0.isconj()) m0.ConjugateSelf();
        } else if (m0.stepj() == 1 && m0.stepi() == N) {
          DoCopy<M,N,ColMajor,RowMajor>(m2.cptr(),m0.ptr());
          if (m2.isconj()) m0.ConjugateSelf();
          if (x2 != T(1)) MultXV<M*N>(x2,m0.ptr());
          if (x1 == T(1))
            AddMM_1<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
          else 
            AddMM<M,N,S1,RowMajor>(x1,m1.cptr(),m0.ptr());
          if (m0.isconj()) m0.ConjugateSelf();
        } else {
          AddMM(x1,m1.View(),x2,m2.View(),m0); 
        }
      } else if (m2.stepj() == 1 && m2.stepi() == N &&
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        if (m0.stepi() == 1 && m0.stepj() == M) {
          DoCopy<M,N,RowMajor,ColMajor>(m2.cptr(),m0.ptr());
          if (m2.isconj()) m0.ConjugateSelf();
          if (x2 != T(1)) MultXV<M*N>(x2,m0.ptr());
          if (x1 == T(1))
            AddMM_1<M,N,S1,ColMajor>(m1.cptr(),m0.ptr());
          else 
            AddMM<M,N,S1,ColMajor>(x1,m1.cptr(),m0.ptr());
          if (m0.isconj()) m0.ConjugateSelf();
        } else if (m0.stepj() == 1 && m0.stepi() == N) {
          DoCopy<M,N,RowMajor,RowMajor>(m2.cptr(),m0.ptr());
          if (m2.isconj()) m0.ConjugateSelf();
          if (x2 != T(1)) MultXV<M*N>(x2,m0.ptr());
          if (x1 == T(1))
            AddMM_1<M,N,S1,RowMajor>(m1.cptr(),m0.ptr());
          else 
            AddMM<M,N,S1,RowMajor>(x1,m1.cptr(),m0.ptr());
          if (m0.isconj()) m0.ConjugateSelf();
        } else {
          AddMM(x1,m1.View(),x2,m2.View(),m0); 
        }
      } else {
        AddMM(x1,m1.View(),x2,m2.View(),m0); 
      }
    }
  private:
    const T x1;
    const SmallMatrix<T1,M,N,S1,I1>& m1;
    const T x2;
    const GenMatrix<T2>& m2;
  };

  template <class T, class T1, class T2, int M, int N, StorageType S2, IndexStyle I2> 
  class SumMm : 
    public MatrixComposite<T> 
  {
  public:
    inline SumMm(const T _x1, const GenMatrix<T1>& _m1, 
        const T _x2, const SmallMatrix<T2,M,N,S2,I2>& _m2) :
      x1(_x1), m1(_m1), x2(_x2), m2(_m2) {}
    inline size_t colsize() const { return M; }
    inline size_t rowsize() const { return N; }
    inline StorageType stor() const { return S2; }
    inline T GetX1() const { return x1; }
    inline const GenMatrix<T1>& GetM1() const { return m1; }
    inline T GetX2() const { return x2; }
    inline const SmallMatrix<T2,M,N,S2,I2>& GetM2() const { return m2; }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    {
      TMVAssert(IsReal(T()));
      if (m1.stepi() == 1 && m1.stepj() == M && 
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        if (m0.stepi() == 1 && m0.stepj() == M) {
          DoCopy<M,N,ColMajor,ColMajor>(m1.cptr(),m0.ptr());
          if (x1 != T(1)) MultXV<M*N>(x1,m0.ptr());
          if (x2 == T(1))
            AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
          else 
            AddMM<M,N,S2,ColMajor>(x2,m2.cptr(),m0.ptr());
        } else if (m0.stepj() == 1 && m0.stepi() == N) {
          DoCopy<M,N,ColMajor,RowMajor>(m1.cptr(),m0.ptr());
          if (x1 != T(1)) MultXV<M*N>(x1,m0.ptr());
          if (x2 == T(1))
            AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
          else 
            AddMM<M,N,S2,RowMajor>(x2,m2.cptr(),m0.ptr());
        } else {
          AddMM(x1,m1.View(),x2,m2.View(),m0); 
        }
      } else if (m2.stepj() == 1 && m2.stepi() == N && 
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        if (m0.stepi() == 1 && m0.stepj() == M) {
          DoCopy<M,N,RowMajor,ColMajor>(m1.cptr(),m0.ptr());
          if (x1 != T(1)) MultXV<M*N>(x1,m0.ptr());
          if (x2 == T(1))
            AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
          else 
            AddMM<M,N,S2,ColMajor>(x2,m2.cptr(),m0.ptr());
        } else if (m0.stepj() == 1 && m0.stepi() == N) {
          DoCopy<M,N,RowMajor,RowMajor>(m1.cptr(),m0.ptr());
          if (x1 != T(1)) MultXV<M*N>(x1,m0.ptr());
          if (x2 == T(1))
            AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
          else 
            AddMM<M,N,S2,RowMajor>(x2,m2.cptr(),m0.ptr());
        } else {
          AddMM(x1,m1.View(),x2,m2.View(),m0); 
        }
      } else {
        AddMM(x1,m1.View(),x2,m2.View(),m0); 
      }
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    {
      if (m1.stepi() == 1 && m1.stepj() == M && 
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        if (m0.stepi() == 1 && m0.stepj() == M) {
          DoCopy<M,N,ColMajor,ColMajor>(m1.cptr(),m0.ptr());
          if (m1.isconj()) m0.ConjugateSelf();
          if (x1 != T(1)) MultXV<M*N>(x1,m0.ptr());
          if (x2 == T(1))
            AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
          else 
            AddMM<M,N,S2,ColMajor>(x2,m2.cptr(),m0.ptr());
          if (m0.isconj()) m0.ConjugateSelf();
        } else if (m0.stepj() == 1 && m0.stepi() == N) {
          DoCopy<M,N,ColMajor,RowMajor>(m1.cptr(),m0.ptr());
          if (m1.isconj()) m0.ConjugateSelf();
          if (x1 != T(1)) MultXV<M*N>(x1,m0.ptr());
          if (x2 == T(1))
            AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
          else 
            AddMM<M,N,S2,RowMajor>(x2,m2.cptr(),m0.ptr());
          if (m0.isconj()) m0.ConjugateSelf();
        } else {
          AddMM(x1,m1.View(),x2,m2.View(),m0); 
        }
      } else if (m2.stepj() == 1 && m2.stepi() == N && 
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        if (m0.stepi() == 1 && m0.stepj() == M) {
          DoCopy<M,N,RowMajor,ColMajor>(m1.cptr(),m0.ptr());
          if (m1.isconj()) m0.ConjugateSelf();
          if (x1 != T(1)) MultXV<M*N>(x1,m0.ptr());
          if (x2 == T(1))
            AddMM_1<M,N,S2,ColMajor>(m2.cptr(),m0.ptr());
          else 
            AddMM<M,N,S2,ColMajor>(x2,m2.cptr(),m0.ptr());
          if (m0.isconj()) m0.ConjugateSelf();
        } else if (m0.stepj() == 1 && m0.stepi() == N) {
          DoCopy<M,N,RowMajor,RowMajor>(m1.cptr(),m0.ptr());
          if (m1.isconj()) m0.ConjugateSelf();
          if (x1 != T(1)) MultXV<M*N>(x1,m0.ptr());
          if (x2 == T(1))
            AddMM_1<M,N,S2,RowMajor>(m2.cptr(),m0.ptr());
          else 
            AddMM<M,N,S2,RowMajor>(x2,m2.cptr(),m0.ptr());
          if (m0.isconj()) m0.ConjugateSelf();
        } else {
          AddMM(x1,m1.View(),x2,m2.View(),m0); 
        }
      } else {
        AddMM(x1,m1.View(),x2,m2.View(),m0); 
      }
    }
  private:
    const T x1;
    const GenMatrix<T1>& m1;
    const T x2;
    const SmallMatrix<T2,M,N,S2,I2>& m2;
  };

  // m += m
  template <class T, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2>
  inline SmallMatrix<T,M,N,S1,I1>& operator+=(
      SmallMatrix<T,M,N,S1,I1>& m1, const SmallMatrix<T,M,N,S2,I2>& m2) 
  { AddMM_1<M,N,S2,S1>(m2.cptr(),m1.ptr()); return m1; }

  template <class T, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2>
  inline SmallMatrix<CT,M,N,S1,I1>& operator+=(
      SmallMatrix<CT,M,N,S1,I1>& m1, const SmallMatrix<T,M,N,S2,I2>& m2) 
  { AddMM_1<M,N,S2,S1>(m2.cptr(),m1.ptr()); return m1; }

  template <class T, int M, int N, StorageType S2, IndexStyle I2>
  inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m1, const SmallMatrix<T,M,N,S2,I2>& m2) 
  {
    TMVAssert(m1.colsize()==M);
    TMVAssert(m1.rowsize()==N);
    AddMM(T(1),m2.View(),m1); 
    return m1; 
  }

  template <class T, int M, int N, StorageType S2, IndexStyle I2>
  inline const MatrixView<CT>& operator+=(
      const MatrixView<CT>& m1, const SmallMatrix<T,M,N,S2,I2>& m2) 
  {
    TMVAssert(m1.colsize()==M);
    TMVAssert(m1.rowsize()==N);
    AddMM(CT(1),m2.View(),m1); 
    return m1; 
  }

  template <class T, int M, int N, StorageType S1, IndexStyle I1>
  inline SmallMatrix<T,M,N,S1,I1>& operator+=(
      SmallMatrix<T,M,N,S1,I1>& m1, const GenMatrix<T>& m2) 
  {
    TMVAssert(m2.colsize()==M);
    TMVAssert(m2.rowsize()==N);
    AddMM(T(1),m2,m1.View()); 
    return m1; 
  }

  template <class T, int M, int N, StorageType S1, IndexStyle I1>
  inline SmallMatrix<CT,M,N,S1,I1>& operator+=(
      SmallMatrix<CT,M,N,S1,I1>& m1, const GenMatrix<T>& m2) 
  {
    TMVAssert(m2.colsize()==M);
    TMVAssert(m2.rowsize()==N);
    AddMM(CT(1),m2,m1.View()); 
    return m1; 
  }

  // m -= m
  template <class T, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2>
  inline SmallMatrix<T,M,N,S1,I1>& operator-=(
      SmallMatrix<T,M,N,S1,I1>& m1, const SmallMatrix<T,M,N,S2,I2>& m2) 
  { AddMM_m1<M,N,S2,S1>(m2.cptr(),m1.ptr()); return m1; }

  template <class T, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2>
  inline SmallMatrix<CT,M,N,S1,I1>& operator-=(
      SmallMatrix<CT,M,N,S1,I1>& m1, const SmallMatrix<T,M,N,S2,I2>& m2) 
  { AddMM_m1<M,N,S2,S1>(m2.cptr(),m1.ptr()); return m1; }

  template <class T, int M, int N, StorageType S2, IndexStyle I2>
  inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m1, const SmallMatrix<T,M,N,S2,I2>& m2) 
  { 
    TMVAssert(m1.colsize()==M);
    TMVAssert(m1.rowsize()==N);
    AddMM(T(-1),m2.View(),m1); 
    return m1; 
  }

  template <class T, int M, int N, StorageType S2, IndexStyle I2>
  inline const MatrixView<CT>& operator-=(
      const MatrixView<CT>& m1, const SmallMatrix<T,M,N,S2,I2>& m2) 
  {
    TMVAssert(m1.colsize()==M);
    TMVAssert(m1.rowsize()==N);
    AddMM(CT(-1),m2.View(),m1); 
    return m1; 
  }

  template <class T, int M, int N, StorageType S1, IndexStyle I1>
  inline SmallMatrix<T,M,N,S1,I1>& operator-=(
      SmallMatrix<T,M,N,S1,I1>& m1, const GenMatrix<T>& m2) 
  {
    TMVAssert(m2.colsize()==M);
    TMVAssert(m2.rowsize()==N);
    AddMM(T(-1),m2,m1.View()); 
    return m1; 
  }

  template <class T, int M, int N, StorageType S1, IndexStyle I1>
  inline SmallMatrix<CT,M,N,S1,I1>& operator-=(
      SmallMatrix<CT,M,N,S1,I1>& m1, const GenMatrix<T>& m2) 
  {
    TMVAssert(m2.colsize()==M);
    TMVAssert(m2.rowsize()==N);
    AddMM(CT(-1),m2,m1.View()); 
    return m1; 
  }

  // m += x*m
  template <class T, class T2, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2>
  inline SmallMatrix<T,M,N,S1,I1>& operator+=(
      SmallMatrix<T,M,N,S1,I1>& m, const ProdXm<T,T2,M,N,S2,I2>& pxm)
  { AddMM<M,N,S2,S1>(pxm.GetX(),pxm.GetM().cptr(),m.ptr()); return m; }

  template <class T, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2>
  inline SmallMatrix<CT,M,N,S1,I1>& operator+=(
      SmallMatrix<CT,M,N,S1,I1>& m, const ProdXm<T,T,M,N,S2,I2>& pxm)
  { AddMM<M,N,S2,S1>(pxm.GetX(),pxm.GetM().cptr(),m.ptr()); return m; }

  template <class T, class T2, int M, int N, StorageType S2, IndexStyle I2>
  inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m, const ProdXm<T,T2,M,N,S2,I2>& pxm)
  {
    TMVAssert(m.colsize()==M);
    TMVAssert(m.rowsize()==N);
    AddMM(pxm.GetX(),pxm.GetM().View(),m); 
    return m; 
  }

  template <class T, int M, int N, StorageType S2, IndexStyle I2>
  inline const MatrixView<CT>& operator+=(
      const MatrixView<CT>& m, const ProdXm<T,T,M,N,S2,I2>& pxm)
  {
    TMVAssert(m.colsize()==M);
    TMVAssert(m.rowsize()==N);
    AddMM(pxm.GetX(),pxm.GetM().View(),m); 
    return m; 
  }

  template <class T, class T2, int M, int N, StorageType S1, IndexStyle I1>
  inline SmallMatrix<T,M,N,S1,I1>& operator+=(
      SmallMatrix<T,M,N,S1,I1>& m, const ProdXM<T,T2>& pxm)
  {
    TMVAssert(pxm.colsize()==M);
    TMVAssert(pxm.rowsize()==N);
    AddMM(pxm.GetX(),pxm.GetM(),m.View()); 
    return m; 
  }

  template <class T, int M, int N, StorageType S1, IndexStyle I1>
  inline SmallMatrix<CT,M,N,S1,I1>& operator+=(
      SmallMatrix<CT,M,N,S1,I1>& m, const ProdXM<T,T>& pxm)
  { 
    TMVAssert(pxm.colsize()==M);
    TMVAssert(pxm.rowsize()==N);
    AddMM(pxm.GetX(),pxm.GetM(),m.View()); 
    return m; 
  }

  // m -= xm
  template <class T, class T2, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2>
  inline SmallMatrix<T,M,N,S1,I1>& operator-=(
      SmallMatrix<T,M,N,S1,I1>& m, const ProdXm<T,T2,M,N,S2,I2>& pxm)
  { AddMM<M,N,S2,S1>(-pxm.GetX(),pxm.GetM().cptr(),m.ptr()); return m; }

  template <class T, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2>
  inline SmallMatrix<CT,M,N,S1,I1>& operator-=(
      SmallMatrix<CT,M,N,S1,I1>& m, const ProdXm<T,T,M,N,S2,I2>& pxm)
  { AddMM<M,N,S2,S1>(-pxm.GetX(),pxm.GetM().cptr(),m.ptr()); return m; }

  template <class T, class T2, int M, int N, StorageType S2, IndexStyle I2>
  inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m, const ProdXm<T,T2,M,N,S2,I2>& pxm)
  {
    TMVAssert(m.colsize()==M);
    TMVAssert(m.rowsize()==N);
    AddMM(-pxm.GetX(),pxm.GetM().View(),m); 
    return m; 
  }

  template <class T, int M, int N, StorageType S2, IndexStyle I2>
  inline const MatrixView<CT>& operator-=(
      const MatrixView<CT>& m, const ProdXm<T,T,M,N,S2,I2>& pxm)
  {
    TMVAssert(m.colsize()==M);
    TMVAssert(m.rowsize()==N);
    AddMM(-pxm.GetX(),pxm.GetM().View(),m); 
    return m; 
  }

  template <class T, class T2, int M, int N, StorageType S1, IndexStyle I1>
  inline SmallMatrix<T,M,N,S1,I1>& operator-=(
      SmallMatrix<T,M,N,S1,I1>& m, const ProdXM<T,T2>& pxm)
  { 
    TMVAssert(pxm.colsize()==M);
    TMVAssert(pxm.rowsize()==N);
    AddMM(-pxm.GetX(),pxm.GetM(),m.View()); 
    return m; 
  }

  template <class T, int M, int N, StorageType S1, IndexStyle I1>
  inline SmallMatrix<CT,M,N,S1,I1>& operator-=(
      SmallMatrix<CT,M,N,S1,I1>& m, const ProdXM<T,T>& pxm)
  { 
    TMVAssert(pxm.colsize()==M);
    TMVAssert(pxm.rowsize()==N);
    AddMM(-pxm.GetX(),pxm.GetM(),m.View()); 
    return m; 
  }


#define GENMATRIX1 SmallMatrix
#define GENMATRIX2 SmallMatrix
#define PRODXM1 ProdXm
#define PRODXM2 ProdXm
#define SUMMM Summm
#define SUMMM_1_1 Summm_1_1
#define SUMMM_1_m1 Summm_1_m1
#define SUMMM_1_x Summm_1_x
#define SUMMM_x_1 Summm_x_1
#define SUMMM_x_m1 Summm_x_m1
#define X1 ,M,N,S1,I1
#define X2 ,M,N,S2,I2
#define X3 ,M,N,S1,S2,I1,I2
#define Y ,int M,int N,StorageType S1,StorageType S2,IndexStyle I1,IndexStyle I2
#include "tmv/TMV_AuxSumMM.h"
#define SUMMM_1_1 Summm_1_1
#define SUMMM_1_m1 Summm_1_m1
#define SUMMM_1_x Summm_1_x
#define SUMMM_x_1 Summm_x_1
#define SUMMM_x_m1 Summm_x_m1
#define X3 ,M,N,S1,S2,I1,I2
#define Y ,int M,int N,StorageType S1,StorageType S2,IndexStyle I1,IndexStyle I2
#include "tmv/TMV_AuxSumMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef SUMMM
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 SmallMatrix
#define GENMATRIX2 GenMatrix
#define SUMMM SummM
#define PRODXM1 ProdXm
#define PRODXM2 ProdXM
#define X1 ,M,N,S1,I1
#define X3 ,M,N,S1,I1
#define Y ,int M,int N, StorageType S1, IndexStyle I1
#include "tmv/TMV_AuxSumMM.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef SUMMM
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 SmallMatrix
#define SUMMM SumMm
#define PRODXM1 ProdXM
#define PRODXM2 ProdXm
#define X2 ,M,N,S2,I2
#define X3 ,M,N,S2,I2
#define Y ,int M,int N, StorageType S2, IndexStyle I2
#include "tmv/TMV_AuxSumMM.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef SUMMM
#undef PRODXM1
#undef PRODXM2


  //
  // Matrix * Matrix
  //

  template <class T, class T1, class T2, int M, int N, int K, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
  class Prodmm_1 : 
    public SmallMatrixComposite<T,M,N>
  {
  public:
    inline Prodmm_1(
        const T DEBUGPARAM(_x), const SmallMatrix<T1,M,K,S1,I1>& _m1, 
        const SmallMatrix<T2,K,N,S2,I2>& _m2) :
      m1(_m1), m2(_m2) { TMVAssert(_x == T(1)); }
    inline StorageType stor() const 
    { return S2==RowMajor ? RowMajor : ColMajor; }
    inline const SmallMatrix<T1,M,K,S1,I1>& GetM1() const { return m1; }
    inline const SmallMatrix<T2,K,N,S2,I2>& GetM2() const { return m2; }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,CStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      MultMM_1<M,N,K,S1,S2,RowMajor>(m1.cptr(), m2.cptr(), m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,CStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      MultMM_1<M,N,K,S1,S2,ColMajor>(m1.cptr(), m2.cptr(), m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,CStyle>& m0) const
    { MultMM_1<M,N,K,S1,S2,RowMajor>(m1.cptr(), m2.cptr(), m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,CStyle>& m0) const
    { MultMM_1<M,N,K,S1,S2,ColMajor>(m1.cptr(), m2.cptr(), m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,FortranStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      MultMM_1<M,N,K,S1,S2,RowMajor>(m1.cptr(), m2.cptr(), m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,FortranStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      MultMM_1<M,N,K,S1,S2,ColMajor>(m1.cptr(), m2.cptr(), m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,FortranStyle>& m0) const
    { MultMM_1<M,N,K,S1,S2,RowMajor>(m1.cptr(), m2.cptr(), m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,FortranStyle>& m0) const
    { MultMM_1<M,N,K,S1,S2,ColMajor>(m1.cptr(), m2.cptr(), m0.ptr()); }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    {
      TMVAssert(IsReal(T()));
      if (m0.stepi() == 1 && m0.stepj() == M && 
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        MultMM_1<M,N,K,S1,S2,ColMajor>(m1.cptr(), m2.cptr(), m0.ptr());
      } else if (m0.stepj() == 1 && m0.stepi() == N &&
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        MultMM_1<M,N,K,S1,S2,RowMajor>(m1.cptr(), m2.cptr(), m0.ptr());
      } else {
        MultMM<false>(T(1), m1.View(), m2.View(), m0); 
      }
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    {
      if (m0.stepi() == 1 && m0.stepj() == M && 
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        MultMM_1<M,N,K,S1,S2,ColMajor>(m1.cptr(), m2.cptr(), m0.ptr());
        if (m0.isconj()) m0.ConjugateSelf();
      } else if (m0.stepj() == 1 && m0.stepi() == N &&
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        MultMM_1<M,N,K,S1,S2,RowMajor>(m1.cptr(), m2.cptr(), m0.ptr());
        if (m0.isconj()) m0.ConjugateSelf();
      } else {
        MultMM<false>(T(1), m1.View(), m2.View(), m0); 
      }
    }
  private:
    const SmallMatrix<T1,M,K,S1,I1>& m1;
    const SmallMatrix<T2,K,N,S2,I2>& m2;
  };

  template <class T, class T1, class T2, int M, int N, int K, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
  class Prodmm : 
    public SmallMatrixComposite<T,M,N>
  {
  public:
    inline Prodmm(const T _x, const SmallMatrix<T1,M,K,S1,I1>& _m1, 
        const SmallMatrix<T2,K,N,S2,I2>& _m2) :
      x(_x), m1(_m1), m2(_m2) {}
    inline StorageType stor() const 
    { return S2==RowMajor ? RowMajor : ColMajor; }
    inline T GetX() const { return x; }
    inline const SmallMatrix<T1,M,K,S1,I1>& GetM1() const { return m1; }
    inline const SmallMatrix<T2,K,N,S2,I2>& GetM2() const { return m2; }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,CStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      MultMM<M,N,K,S1,S2,RowMajor>(x,m1.cptr(), m2.cptr(), m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,CStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      MultMM<M,N,K,S1,S2,ColMajor>(x,m1.cptr(), m2.cptr(), m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,CStyle>& m0) const
    { MultMM<M,N,K,S1,S2,RowMajor>(x,m1.cptr(), m2.cptr(), m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,CStyle>& m0) const
    { MultMM<M,N,K,S1,S2,ColMajor>(x,m1.cptr(), m2.cptr(), m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,FortranStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      MultMM<M,N,K,S1,S2,RowMajor>(x,m1.cptr(), m2.cptr(), m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,FortranStyle>& m0) const
    {
      TMVAssert(IsReal(T()));
      MultMM<M,N,K,S1,S2,ColMajor>(x,m1.cptr(), m2.cptr(), m0.ptr());
    }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,FortranStyle>& m0) const
    { MultMM<M,N,K,S1,S2,RowMajor>(x,m1.cptr(), m2.cptr(), m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,FortranStyle>& m0) const
    { MultMM<M,N,K,S1,S2,ColMajor>(x,m1.cptr(), m2.cptr(), m0.ptr()); }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    {
      TMVAssert(IsReal(T()));
      if (m0.stepi() == 1 && m0.stepj() == M && 
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        MultMM<M,N,K,S1,S2,ColMajor>(x,m1.cptr(), m2.cptr(), m0.ptr());
      } else if (m0.stepj() == 1 && m0.stepi() == N &&
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        MultMM<M,N,K,S1,S2,RowMajor>(x,m1.cptr(), m2.cptr(), m0.ptr());
      } else {
        MultMM<false>(x, m1.View(), m2.View(), m0); 
      }
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    {
      if (m0.stepi() == 1 && m0.stepj() == M && 
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        MultMM<M,N,K,S1,S2,ColMajor>(x,m1.cptr(), m2.cptr(), m0.ptr());
        if (m0.isconj()) m0.ConjugateSelf();
      } else if (m0.stepj() == 1 && m0.stepi() == N &&
          !SameStorage(m0,m1) && !SameStorage(m0,m2)) {
        MultMM<M,N,K,S1,S2,RowMajor>(x,m1.cptr(), m2.cptr(), m0.ptr());
        if (m0.isconj()) m0.ConjugateSelf();
      } else {
        MultMM<false>(x, m1.View(), m2.View(), m0); 
      }
    }
  private:
    const T x;
    const SmallMatrix<T1,M,K,S1,I1>& m1;
    const SmallMatrix<T2,K,N,S2,I2>& m2;
  };

  template <class T, class T1, class T2, int K, int N, StorageType S2, IndexStyle I2> 
  class ProdMm : 
    public MatrixComposite<T>
  {
  public:
    inline ProdMm(const T _x, const GenMatrix<T1>& _m1, 
        const SmallMatrix<T2,K,N,S2,I2>& _m2) :
      x(_x), m1(_m1), m2(_m2) { TMVAssert(m1.rowsize() == K); }
    inline size_t colsize() const { return m1.colsize(); }
    inline size_t rowsize() const { return N; }
    inline StorageType stor() const 
    { return S2==RowMajor ? RowMajor : ColMajor; }
    inline T GetX() const { return x; }
    inline const GenMatrix<T1>& GetM1() const { return m1; }
    inline const SmallMatrix<T2,K,N,S2,I2>& GetM2() const { return m2; }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    { TMVAssert(IsReal(T())); MultMM<false>(x, m1, m2.View(), m0); }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    { MultMM<false>(x, m1, m2.View(), m0); }
  private:
    const T x;
    const GenMatrix<T1>& m1;
    const SmallMatrix<T2,K,N,S2,I2>& m2;
  };

  template <class T, class T1, class T2, int M, int K, StorageType S1, IndexStyle I1> 
  class ProdmM : 
    public MatrixComposite<T>
  {
  public:
    inline ProdmM(const T _x, const SmallMatrix<T1,M,K,S1,I1>& _m1, 
        const GenMatrix<T2>& _m2) :
      x(_x), m1(_m1), m2(_m2) { TMVAssert(m2.colsize() == K); }
    inline size_t colsize() const { return M; }
    inline size_t rowsize() const { return m2.rowsize(); }
    inline StorageType stor() const 
    { return m2.stor()==RowMajor ? RowMajor : ColMajor; }
    inline T GetX() const { return x; }
    inline const SmallMatrix<T1,M,K,S1,I1>& GetM1() const { return m1; }
    inline const GenMatrix<T2>& GetM2() const { return m2; }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    { TMVAssert(IsReal(T())); MultMM<false>(x, m1.View(), m2, m0); }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    { MultMM<false>(x, m1.View(), m2, m0); }
  private:
    const T x;
    const SmallMatrix<T1,M,K,S1,I1>& m1;
    const GenMatrix<T2>& m2;
  };

  // m *= m
  template <class T, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2>
  inline SmallMatrix<T,M,N,S1,I1>& operator*=(
      SmallMatrix<T,M,N,S1,I1>& m1, const SmallMatrix<T,N,N,S2,I2>& m2)
  { 
    SmallMatrix<T,M,N,S1> m1_copy(m1);
    MultMM_1<M,N,N,S1,S2,S1>(m1_copy.cptr(),m2.cptr(),m1.ptr());
    return m1;
  }

  template <class T, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2>
  inline SmallMatrix<CT,M,N,S1,I1>& operator*=(
      SmallMatrix<CT,M,N,S1,I1>& m1, const SmallMatrix<T,N,N,S2,I2>& m2)
  {
    SmallMatrix<CT,M,N,S1> m1_copy(m1);
    MultMM_1<M,N,N,S1,S2,S1>(m1_copy.cptr(),m2.cptr(),m1.ptr());
    return m1;
  }

  template <class T, int N, StorageType S2, IndexStyle I2>
  inline const MatrixView<T>& operator*=(
      const MatrixView<T>& m1, const SmallMatrix<T,N,N,S2,I2>& m2)
  { 
    TMVAssert(m1.rowsize() == N);
    MultMM<false>(T(1),m1,m2.View(),m1); 
    return m1; 
  }

  template <class T, int N, StorageType S2, IndexStyle I2>
  inline const MatrixView<CT>& operator*=(
      const MatrixView<CT>& m1, const SmallMatrix<T,N,N,S2,I2>& m2)
  { 
    TMVAssert(m1.rowsize() == N);
    MultMM<false>(T(1),m1,m2.View(),m1); 
    return m1; 
  }

  template <class T, int M, int N, StorageType S1, IndexStyle I1>
  inline SmallMatrix<T,M,N,S1,I1>& operator*=(
      SmallMatrix<T,M,N,S1,I1>& m1, const GenMatrix<T>& m2)
  { 
    TMVAssert(m2.rowsize() == N);
    MultMM<false>(T(1),m1.View(),m2,m1.View()); 
    return m1; 
  }

  template <class T, int M, int N, StorageType S1, IndexStyle I1>
  inline SmallMatrix<CT,M,N,S1,I1>& operator*=(
      SmallMatrix<CT,M,N,S1,I1>& m1, const GenMatrix<T>& m2)
  {
    TMVAssert(m2.rowsize() == N);
    MultMM<false>(T(1),m1.View(),m2,m1.View()); 
    return m1; 
  }

  // m *= xm
  template <class T, class T2, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2>
  inline SmallMatrix<T,M,N,S1,I1>& operator*=(
      SmallMatrix<T,M,N,S1,I1>& m1, const ProdXm<T,T2,N,N,S2,I2>& pxm)
  { 
    SmallMatrix<T,M,N,S1> m1_copy(m1);
    MultMM<M,N,N,S1,S2,S1>(pxm.GetX(),m1_copy.cptr(),pxm.GetM().cptr(),
        m1.ptr());
    return m1;
  }

  template <class T, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2>
  inline SmallMatrix<CT,M,N,S1,I1>& operator*=(
      SmallMatrix<CT,M,N,S1,I1>& m1, const ProdXm<T,T,N,N,S2,I2>& pxm)
  {
    SmallMatrix<CT,M,N,S1> m1_copy(m1);
    MultMM<M,N,N,S1,S2,S1>(pxm.GetX(),m1_copy.cptr(),pxm.GetM().cptr(),
        m1.ptr());
    return m1;
  }

  template <class T, class T2, int N, StorageType S2, IndexStyle I2>
  inline const MatrixView<T>& operator*=(
      const MatrixView<T>& m1, const ProdXm<T,T2,N,N,S2,I2>& pxm)
  {
    TMVAssert(m1.rowsize() == N);
    MultMM<false>(pxm.GetX(),m1,pxm.GetM().View(),m1);
    return m1; 
  }

  template <class T, int N, StorageType S2, IndexStyle I2>
  inline const MatrixView<CT>& operator*=(
      const MatrixView<CT>& m1, const ProdXm<T,T,N,N,S2,I2>& pxm)
  {
    TMVAssert(m1.rowsize() == N);
    MultMM<false>(pxm.GetX(),m1,pxm.GetM().View(),m1);
    return m1; 
  }

  template <class T, class T2, int M, int N, StorageType S1, IndexStyle I1>
  inline SmallMatrix<T,M,N,S1,I1>& operator*=(
      SmallMatrix<T,M,N,S1,I1>& m1, const ProdXM<T,T2>& pxm)
  { 
    TMVAssert(pxm.colsize() == N);
    TMVAssert(pxm.rowsize() == N);
    MultMM<false>(pxm.GetX(),m1.View(),pxm.GetM(),m1.View()); 
    return m1; 
  }

  template <class T, int M, int N, StorageType S1, IndexStyle I1>
  inline SmallMatrix<CT,M,N,S1,I1>& operator*=(
      SmallMatrix<CT,M,N,S1,I1>& m1, const ProdXM<T,T>& pxm)
  { 
    TMVAssert(pxm.colsize() == N);
    TMVAssert(pxm.rowsize() == N);
    MultMM<false>(pxm.GetX(),m1.View(),pxm.GetM(),m1.View()); 
    return m1; 
  }

  // m += mm
  template <class T, class T2, class T3, int M, int N, int K, StorageType S1, StorageType S2, StorageType S3, IndexStyle I1, IndexStyle I2, IndexStyle I3>
  inline SmallMatrix<T,M,N,S1,I1>& operator+=(
      SmallMatrix<T,M,N,S1,I1>& m, 
      const Prodmm_1<T,T2,T3,M,N,K,S2,S3,I2,I3>& pmm)
  { 
    AddMultMM_1<M,N,K,S2,S3,S1>(pmm.GetM1().cptr(),pmm.GetM2().cptr(),
        m.ptr());
    return m;
  }

  template <class T, int M, int N, int K, StorageType S1, StorageType S2, StorageType S3, IndexStyle I1, IndexStyle I2, IndexStyle I3>
  inline SmallMatrix<CT,M,N,S1,I1>& operator+=(
      SmallMatrix<CT,M,N,S1,I1>& m, 
      const Prodmm_1<T,T,T,M,N,K,S2,S3,I2,I3>& pmm)
  {
    AddMultMM_1<M,N,K,S2,S3,S1>(pmm.GetM1().cptr(),pmm.GetM2().cptr(),
        m.ptr());
    return m;
  }

  template <class T, class T2, class T3, int M, int N, int K, StorageType S2, StorageType S3, IndexStyle I2, IndexStyle I3>
  inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m, const Prodmm_1<T,T2,T3,M,N,K,S2,S3,I2,I3>& pmm)
  { 
    TMVAssert(m.colsize() == M);
    TMVAssert(m.rowsize() == N);
    MultMM<true>(T(1),pmm.GetM1().View(),pmm.GetM2().View(),m);
    return m; 
  }

  template <class T, int M, int N, int K, StorageType S2, StorageType S3, IndexStyle I2, IndexStyle I3>
  inline const MatrixView<CT>& operator+=(
      const MatrixView<CT>& m, const Prodmm_1<T,T,T,M,N,K,S2,S3,I2,I3>& pmm)
  { 
    TMVAssert(m.colsize() == M);
    TMVAssert(m.rowsize() == N);
    MultMM<true>(T(1),pmm.GetM1().View(),pmm.GetM2().View(),m);
    return m; 
  }

  // m -= mm
  template <class T, class T2, class T3, int M, int N, int K, StorageType S1, StorageType S2, StorageType S3, IndexStyle I1, IndexStyle I2, IndexStyle I3>
  inline SmallMatrix<T,M,N,S1,I1>& operator-=(
      SmallMatrix<T,M,N,S1,I1>& m, 
      const Prodmm_1<T,T2,T3,M,N,K,S2,S3,I2,I3>& pmm)
  { 
    AddMultMM_m1<M,N,K,S2,S3,S1>(pmm.GetM1().cptr(),pmm.GetM2().cptr(),
        m.ptr());
    return m; 
  }

  template <class T, int M, int N, int K, StorageType S1, StorageType S2, StorageType S3, IndexStyle I1, IndexStyle I2, IndexStyle I3>
  inline SmallMatrix<CT,M,N,S1,I1>& operator-=(
      SmallMatrix<CT,M,N,S1,I1>& m, 
      const Prodmm_1<T,T,T,M,N,K,S2,S3,I2,I3>& pmm)
  { 
    AddMultMM_m1<M,N,K,S2,S3,S1>(pmm.GetM1().cptr(),pmm.GetM2().cptr(),
        m.ptr());
    return m; 
  }

  template <class T, class T2, class T3, int M, int N, int K, StorageType S2, StorageType S3, IndexStyle I2, IndexStyle I3>
  inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m, const Prodmm_1<T,T2,T3,M,N,K,S2,S3,I2,I3>& pmm)
  { 
    TMVAssert(m.colsize() == M);
    TMVAssert(m.rowsize() == N);
    MultMM<true>(T(-1),pmm.GetM1().View(),pmm.GetM2().View(),m);
    return m; 
  }

  template <class T, int M, int N, int K, StorageType S2, StorageType S3, IndexStyle I2, IndexStyle I3>
  inline const MatrixView<CT>& operator-=(
      const MatrixView<CT>& m, const Prodmm_1<T,T,T,M,N,K,S2,S3,I2,I3>& pmm)
  { 
    TMVAssert(m.colsize() == M);
    TMVAssert(m.rowsize() == N);
    MultMM<true>(T(-1),pmm.GetM1().View(),pmm.GetM2().View(),m);
    return m; 
  }

  // m += xmm
  template <class T, class T2, class T3, int M, int N, int K, StorageType S1, StorageType S2, StorageType S3, IndexStyle I1, IndexStyle I2, IndexStyle I3>
  inline SmallMatrix<T,M,N,S1,I1>& operator+=(
      SmallMatrix<T,M,N,S1,I1>& m, 
      const Prodmm<T,T2,T3,M,N,K,S2,S3,I2,I3>& pmm)
  { 
    AddMultMM<M,N,K,S2,S3,S1>(pmm.GetX(),pmm.GetM1().cptr(),
        pmm.GetM2().cptr(),m.ptr());
    return m; 
  }

  template <class T, int M, int N, int K, StorageType S1, StorageType S2, StorageType S3, IndexStyle I1, IndexStyle I2, IndexStyle I3>
  inline SmallMatrix<CT,M,N,S1,I1>& operator+=(
      SmallMatrix<CT,M,N,S1,I1>& m, 
      const Prodmm<T,T,T,M,N,K,S2,S3,I2,I3>& pmm)
  {
    AddMultMM<M,N,K,S2,S3,S1>(pmm.GetX(),pmm.GetM1().cptr(),
        pmm.GetM2().cptr(),m.ptr());
    return m;
  }

  template <class T, class T2, class T3, int M, int N, int K, StorageType S2, StorageType S3, IndexStyle I2, IndexStyle I3>
  inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m, const Prodmm<T,T2,T3,M,N,K,S2,S3,I2,I3>& pmm)
  {
    TMVAssert(m.colsize() == M);
    TMVAssert(m.rowsize() == N);
    MultMM<true>(pmm.GetX(),pmm.GetM1().View(),pmm.GetM2().View(),m);
    return m; 
  }

  template <class T, int M, int N, int K, StorageType S2, StorageType S3, IndexStyle I2, IndexStyle I3>
  inline const MatrixView<CT>& operator+=(
      const MatrixView<CT>& m, const Prodmm<T,T,T,M,N,K,S2,S3,I2,I3>& pmm)
  {
    TMVAssert(m.colsize() == M);
    TMVAssert(m.rowsize() == N);
    MultMM<true>(pmm.GetX(),pmm.GetM1().View(),pmm.GetM2().View(),m);
    return m; 
  }

  template <class T, class T2, class T3, int M, int N, int K, StorageType S1, IndexStyle I1>
  inline SmallMatrix<T,M,N,S1,I1>& operator+=(
      SmallMatrix<T,M,N,S1,I1>& m, const ProdMM<T,T2,T3>& pmm)
  {
    TMVAssert(pmm.colsize() == M);
    TMVAssert(pmm.rowsize() == N);
    MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m.View()); 
    return m; 
  }

  template <class T, int M, int N, int K, StorageType S1, IndexStyle I1>
  inline SmallMatrix<CT,M,N,S1,I1>& operator+=(
      SmallMatrix<CT,M,N,S1,I1>& m, const ProdMM<T,T,T>& pmm)
  {
    TMVAssert(pmm.colsize() == M);
    TMVAssert(pmm.rowsize() == N);
    MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m.View()); 
    return m; 
  }

  // m -= xmm
  template <class T, class T2, class T3, int M, int N, int K, StorageType S1, StorageType S2, StorageType S3, IndexStyle I1, IndexStyle I2, IndexStyle I3>
  inline SmallMatrix<T,M,N,S1,I1>& operator-=(
      SmallMatrix<T,M,N,S1,I1>& m, 
      const Prodmm<T,T2,T3,M,N,K,S2,S3,I2,I3>& pmm)
  { 
    AddMultMM<M,N,K,S2,S3,S1>(-pmm.GetX(),pmm.GetM1().cptr(),
        pmm.GetM2().cptr(),m.ptr());
    return m;
  }

  template <class T, int M, int N, int K, StorageType S1, StorageType S2, StorageType S3, IndexStyle I1, IndexStyle I2, IndexStyle I3>
  inline SmallMatrix<CT,M,N,S1,I1>& operator-=(
      SmallMatrix<CT,M,N,S1,I1>& m, 
      const Prodmm<T,T,T,M,N,K,S2,S3,I2,I3>& pmm)
  { 
    AddMultMM<M,N,K,S2,S3,S1>(-pmm.GetX(),pmm.GetM1().cptr(),
        pmm.GetM2().cptr(),m.ptr());
    return m;
  }

  template <class T, class T2, class T3, int M, int N, int K, StorageType S2, StorageType S3, IndexStyle I2, IndexStyle I3>
  inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m, const Prodmm<T,T2,T3,M,N,K,S2,S3,I2,I3>& pmm)
  {
    TMVAssert(m.colsize() == M);
    TMVAssert(m.rowsize() == N);
    MultMM<true>(-pmm.GetX(),pmm.GetM1().View(),pmm.GetM2().View(),m);
    return m; 
  }

  template <class T, int M, int N, int K, StorageType S2, StorageType S3, IndexStyle I2, IndexStyle I3>
  inline const MatrixView<CT>& operator-=(
      const MatrixView<CT>& m, const Prodmm<T,T,T,M,N,K,S2,S3,I2,I3>& pmm)
  {
    TMVAssert(m.colsize() == M);
    TMVAssert(m.rowsize() == N);
    MultMM<true>(-pmm.GetX(),pmm.GetM1().View(),pmm.GetM2().View(),m);
    return m; 
  }

  template <class T, class T2, class T3, int M, int N, int K, StorageType S1, IndexStyle I1>
  inline SmallMatrix<T,M,N,S1,I1>& operator-=(
      SmallMatrix<T,M,N,S1,I1>& m, const ProdMM<T,T2,T3>& pmm)
  { 
    TMVAssert(pmm.colsize() == M);
    TMVAssert(pmm.rowsize() == N);
    MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m.View()); 
    return m; 
  }

  template <class T, int M, int N, int K, StorageType S1, IndexStyle I1>
  inline SmallMatrix<CT,M,N,S1,I1>& operator-=(
      SmallMatrix<CT,M,N,S1,I1>& m, const ProdMM<T,T,T>& pmm)
  {
    TMVAssert(pmm.colsize() == M);
    TMVAssert(pmm.rowsize() == N);
    MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m.View()); 
    return m; 
  }

#define GENMATRIX1 SmallMatrix
#define GENMATRIX2 SmallMatrix
#define PRODMM Prodmm
#define PRODMM_1 Prodmm_1
#define PRODXM1 ProdXm
#define PRODXM2 ProdXm
#define X1 ,M,K,S1,I1
#define X2 ,K,N,S2,I2
#define X3 ,M,N,K,S1,S2,I1,I2
#define Y ,int M,int N,int K,StorageType S1,StorageType S2, IndexStyle I1, IndexStyle I2
#include "tmv/TMV_AuxProdMM.h"
#define PRODMM_1 Prodmm_1
#define X3 ,M,N,K,S1,S2,I1,I2
#define Y ,int M,int N,int K,StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 SmallMatrix
#define GENMATRIX2 GenMatrix
#define PRODMM ProdmM
#define PRODXM1 ProdXm
#define PRODXM2 ProdXM
#define X1 ,M,K,S1,I1
#define X3 ,M,K,S1,I1
#define Y ,int M,int K, StorageType S1,IndexStyle I1
#include "tmv/TMV_AuxProdMM.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 SmallMatrix
#define PRODMM ProdMm
#define PRODXM1 ProdXM
#define PRODXM2 ProdXm
#define X2 ,K,N,S2,I2
#define X3 ,K,N,S2,I2
#define Y ,int N,int K,StorageType S2,IndexStyle I2
#include "tmv/TMV_AuxProdMM.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2

  //
  // Matrix * Vector
  //

  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I1, IndexStyle I2> 
  class Prodmv_1 :
    public SmallVectorComposite<T,M>
  {
  public:
    inline Prodmv_1(
        const T DEBUGPARAM(_x), const SmallMatrix<T1,M,N,S,I1>& _m,
        const SmallVector<T2,N,I2>& _v) :
      m(_m), v(_v) { TMVAssert(_x == T(1)); }
    inline const SmallMatrix<T1,M,N,S,I1>& GetM() const { return m; }
    inline const SmallVector<T2,N,I2>& GetV() const { return v; }
    inline void AssignTov(SmallVector<RealType(T),M,CStyle>& v0) const
    {
      TMVAssert(IsReal(T()));
      MultMV_1<M,N,S>(m.cptr(),v.cptr(),v0.ptr());
    }
    inline void AssignTov(SmallVector<ComplexType(T),M,CStyle>& v0) const
    {
      MultMV_1<M,N,S>(m.cptr(),v.cptr(),v0.ptr());
    }
    inline void AssignTov(SmallVector<RealType(T),M,FortranStyle>& v0) const
    {
      TMVAssert(IsReal(T()));
      MultMV_1<M,N,S>(m.cptr(),v.cptr(),v0.ptr());
    }
    inline void AssignTov(SmallVector<ComplexType(T),M,FortranStyle>& v0) const
    { 
      MultMV_1<M,N,S>(m.cptr(),v.cptr(),v0.ptr());
    }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    {
      TMVAssert(IsReal(T()));
      if (v0.step() == 1) 
        MultMV_1<M,N,S>(m.cptr(),v.cptr(),v0.ptr());
      else
        MultMV<false>(T(1),m.View(),v.View(),v0); 
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    { 
      if (v0.step() == 1) {
        MultMV_1<M,N,S>(m.cptr(),v.cptr(),v0.ptr());
        if (v0.isconj()) v0.ConjugateSelf();
      }
      else
        MultMV<false>(T(1),m.View(),v.View(),v0);
    }
  private:
    const SmallMatrix<T1,M,N,S,I1>& m;
    const SmallVector<T2,N,I2>& v;
  };

  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I1, IndexStyle I2> 
  class Prodmv :
    public SmallVectorComposite<T,M>
  {
  public:
    inline Prodmv(const T _x, const SmallMatrix<T1,M,N,S,I1>& _m,
        const SmallVector<T2,N,I2>& _v) :
      x(_x), m(_m), v(_v) {}
    inline T GetX() const { return x; }
    inline const SmallMatrix<T1,M,N,S,I1>& GetM() const { return m; }
    inline const SmallVector<T2,N,I2>& GetV() const { return v; }
    inline void AssignTov(SmallVector<RealType(T),M,CStyle>& v0) const
    {
      TMVAssert(IsReal(T()));
      MultMV<M,N,S>(x,m.cptr(),v.cptr(),v0.ptr());
    }
    inline void AssignTov(SmallVector<ComplexType(T),M,CStyle>& v0) const
    { MultMV<M,N,S>(x,m.cptr(),v.cptr(),v0.ptr()); }
    inline void AssignTov(SmallVector<RealType(T),M,FortranStyle>& v0) const
    {
      TMVAssert(IsReal(T()));
      MultMV<M,N,S>(x,m.cptr(),v.cptr(),v0.ptr());
    }
    inline void AssignTov(SmallVector<ComplexType(T),M,FortranStyle>& v0) const
    { MultMV<M,N,S>(x,m.cptr(),v.cptr(),v0.ptr()); }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    { 
      TMVAssert(IsReal(T()));
      if (v0.step() == 1) 
        MultMV<M,N,S>(x,m.cptr(),v.cptr(),v0.ptr());
      else
        MultMV<false>(x,m.View(),v.View(),v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    { 
      if (v0.step() == 1) {
        MultMV<M,N,S>(x,m.cptr(),v.cptr(),v0.ptr());
        if (v0.isconj()) v0.ConjugateSelf();
      }
      else
        MultMV<false>(x,m.View(),v.View(),v0);
    }
  private:
    const T x;
    const SmallMatrix<T1,M,N,S,I1>& m;
    const SmallVector<T2,N,I2>& v;
  };

  template <class T, class T1, class T2, int N, IndexStyle I> 
  class ProdMv :
    public VectorComposite<T>
  {
  public:
    inline ProdMv(const T _x, const GenMatrix<T1>& _m,
        const SmallVector<T2,N,I>& _v) :
      x(_x), m(_m), v(_v) { TMVAssert(m.rowsize() == N); }
    inline size_t size() const { return m.colsize(); }
    inline T GetX() const { return x; }
    inline const GenMatrix<T1>& GetM() const { return m; }
    inline const SmallVector<T2,N,I>& GetV() const { return v; }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    { 
      TMVAssert(IsReal(T()));
      MultMV<false>(x,m,v.View(),v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    { MultMV<false>(x,m,v.View(),v0); }
  private:
    const T x;
    const GenMatrix<T1>& m;
    const SmallVector<T2,N,I>& v;
  };

  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I> 
  class ProdmV :
    public VectorComposite<T>
  {
  public:
    inline ProdmV(const T _x, const SmallMatrix<T1,M,N,S,I>& _m,
        const GenVector<T2>& _v) :
      x(_x), m(_m), v(_v) { TMVAssert(v.size() == N); }
    inline size_t size() const { return M; }
    inline T GetX() const { return x; }
    inline const SmallMatrix<T1,M,N,S,I>& GetM() const { return m; }
    inline const GenVector<T2>& GetV() const { return v; }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    { 
      TMVAssert(IsReal(T()));
      if (v.step() == 1 && v0.step() == 1 && !SameStorage(v0,v)) 
        if (x == T(1))
          MultMV_1<M,N,S>(m.cptr(),v.cptr(),v0.ptr());
        else
          MultMV<M,N,S>(x,m.cptr(),v.cptr(),v0.ptr());
      else
        MultMV<false>(x,m.View(),v,v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    {
      if (v.step() == 1 && v0.step() == 1 && !SameStorage(v0,v) && 
          !v.isconj()) {
        if (x == T(1))
          MultMV_1<M,N,S>(m.cptr(),v.cptr(),v0.ptr());
        else
          MultMV<M,N,S>(x,m.cptr(),v.cptr(),v0.ptr());
        if (v0.isconj()) v0.ConjugateSelf();
      }
      else
        MultMV<false>(x,m.View(),v,v0);
    }
  private:
    const T x;
    const SmallMatrix<T1,M,N,S,I>& m;
    const GenVector<T2>& v;
  };


  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I1, IndexStyle I2> 
  class Prodvm_1 :
    public SmallVectorComposite<T,N>
  {
  public:
    inline Prodvm_1(
        const T DEBUGPARAM(_x), const SmallVector<T1,M,I1>& _v,
        const SmallMatrix<T2,M,N,S,I2>& _m) :
      v(_v), m(_m) { TMVAssert(_x == T(1)); }
    inline const SmallVector<T1,M,I1>& GetV() const { return v; }
    inline const SmallMatrix<T2,M,N,S,I2>& GetM() const { return m; }
    inline void AssignTov(SmallVector<RealType(T),N,CStyle>& v0) const
    {
      TMVAssert(IsReal(T()));
      MultVM_1<M,N,S>(v.cptr(),m.cptr(),v0.ptr());
    }
    inline void AssignTov(SmallVector<ComplexType(T),N,CStyle>& v0) const
    { MultVM_1<M,N,S>(v.cptr(),m.cptr(),v0.ptr()); }
    inline void AssignTov(SmallVector<RealType(T),N,FortranStyle>& v0) const
    {
      TMVAssert(IsReal(T()));
      MultVM_1<M,N,S>(v.cptr(),m.cptr(),v0.ptr());
    }
    inline void AssignTov(SmallVector<ComplexType(T),N,FortranStyle>& v0) const
    { MultVM_1<M,N,S>(v.cptr(),m.cptr(),v0.ptr()); }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    { 
      TMVAssert(IsReal(T()));
      if (v0.step() == 1) 
        MultVM_1<M,N,S>(v.cptr(),m.cptr(),v0.ptr());
      else
        MultMV<false>(T(1),m.Transpose(),v.View(),v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    {
      if (v0.step() == 1) {
        MultVM_1<M,N,S>(v.cptr(),m.cptr(),v0.ptr());
        if (v0.isconj()) v0.ConjugateSelf();
      }
      else
        MultMV<false>(T(1),m.Transpose(),v.View(),v0);
    }
  private:
    const SmallVector<T1,M,I1>& v;
    const SmallMatrix<T2,M,N,S,I2>& m;
  };

  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I1, IndexStyle I2> 
  class Prodvm :
    public SmallVectorComposite<T,N>
  {
  public:
    inline Prodvm(const T _x, const SmallVector<T1,M,I1>& _v,
        const SmallMatrix<T2,M,N,S,I2>& _m) :
      x(_x), v(_v), m(_m) { }
    inline T GetX() const { return x; }
    inline const SmallVector<T1,M,I1>& GetV() const { return v; }
    inline const SmallMatrix<T2,M,N,S,I2>& GetM() const { return m; }
    inline void AssignTov(SmallVector<RealType(T),N,CStyle>& v0) const
    { 
      TMVAssert(IsReal(T()));
      MultVM<M,N,S>(x,v.cptr(),m.cptr(),v0.ptr());
    }
    inline void AssignTov(SmallVector<ComplexType(T),N,CStyle>& v0) const
    { MultVM<M,N,S>(x,v.cptr(),m.cptr(),v0.ptr()); }
    inline void AssignTov(SmallVector<RealType(T),N,FortranStyle>& v0) const
    {
      TMVAssert(IsReal(T()));
      MultVM<M,N,S>(x,v.cptr(),m.cptr(),v0.ptr());
    }
    inline void AssignTov(SmallVector<ComplexType(T),N,FortranStyle>& v0) const
    { MultVM<M,N,S>(x,v.cptr(),m.cptr(),v0.ptr()); }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    {
      TMVAssert(IsReal(T()));
      if (v0.step() == 1) 
        MultVM<M,N,S>(x,v.cptr(),m.cptr(),v0.ptr());
      else
        MultMV<false>(x,m.Transpose(),v.View(),v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    {
      if (v0.step() == 1) {
        MultVM<M,N,S>(x,v.cptr(),m.cptr(),v0.ptr());
        if (v0.isconj()) v0.ConjugateSelf();
      }
      else
        MultMV<false>(x,m.Transpose(),v.View(),v0);
    }
  private:
    const T x;
    const SmallVector<T1,M,I1>& v;
    const SmallMatrix<T2,M,N,S,I2>& m;
  };

  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I> 
  class ProdVm :
    public VectorComposite<T>
  {
  public:
    inline ProdVm(const T _x, const GenVector<T1>& _v,
        const SmallMatrix<T2,M,N,S,I>& _m) :
      x(_x), v(_v), m(_m) { TMVAssert(v.size() == M); }
    inline size_t size() const { return N; }
    inline T GetX() const { return x; }
    inline const GenVector<T1>& GetV() const { return v; }
    inline const SmallMatrix<T2,M,N,S,I>& GetM() const { return m; }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    { 
      TMVAssert(IsReal(T()));
      if (v0.step() == 1 && v.step() == 1 && !SameStorage(v0,v)) 
        if (x == T(1))
          MultVM_1<M,N,S>(v.cptr(),m.cptr(),v0.ptr());
        else
          MultVM<M,N,S>(x,v.cptr(),m.cptr(),v0.ptr());
      else
        MultMV<false>(x,m.Transpose(),v,v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    {
      if (v0.step() == 1 && v.step() == 1 && !SameStorage(v0,v) &&
          !v.isconj()) {
        if (x == T(1))
          MultVM_1<M,N,S>(v.cptr(),m.cptr(),v0.ptr());
        else
          MultVM<M,N,S>(x,v.cptr(),m.cptr(),v0.ptr());
        if (v0.isconj()) v0.ConjugateSelf();
      }
      else
        MultMV<false>(x,m.Transpose(),v,v0);
    }
  private:
    const T x;
    const GenVector<T1>& v;
    const SmallMatrix<T2,M,N,S,I>& m;
  };

  template <class T, class T1, class T2, int M, IndexStyle I> 
  class ProdvM :
    public VectorComposite<T>
  {
  public:
    inline ProdvM(const T _x, const SmallVector<T1,M,I>& _v,
        const GenMatrix<T2>& _m) :
      x(_x), v(_v), m(_m) { TMVAssert(m.colsize() == M); }
    inline size_t size() const { return m.rowsize(); }
    inline T GetX() const { return x; }
    inline const SmallVector<T1,M,I>& GetV() const { return v; }
    inline const GenMatrix<T2>& GetM() const { return m; }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    {
      TMVAssert(IsReal(T()));
      MultMV<false>(x,m.Transpose(),v.View(),v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    { MultMV<false>(x,m.Transpose(),v.View(),v0); }
  private:
    const T x;
    const SmallVector<T1,M,I>& v;
    const GenMatrix<T2>& m;
  };

  // v *= m
  template <class T, int N, StorageType S, IndexStyle I1, IndexStyle I2> 
  inline SmallVector<T,N,I1>& operator*=(
      SmallVector<T,N,I1>& v, const SmallMatrix<T,N,N,S,I2>& m)
  { 
    SmallVector<T,N> v_copy(v);
    MultVM_1<N,N,S>(v_copy.cptr(),m.cptr(),v.ptr());
    return v; 
  }

  template <class T, int N, StorageType S, IndexStyle I1, IndexStyle I2> 
  inline SmallVector<CT,N,I1>& operator*=(
      SmallVector<CT,N,I1>& v, const SmallMatrix<T,N,N,S,I2>& m)
  {
    SmallVector<CT,N> v_copy(v);
    MultVM_1<N,N,S>(v_copy.cptr(),m.cptr(),v.ptr());
    return v; 
  }

  template <class T, int N, StorageType S, IndexStyle I2> 
  inline const VectorView<T>& operator*=(
      const VectorView<T>& v, const SmallMatrix<T,N,N,S,I2>& m)
  {
    TMVAssert(v.size() == N);
    if (v.step() == 1 && !v.isconj()) {
      SmallVector<T,N> v_copy(v);
      MultVM_1<N,N,S>(v_copy.cptr(),m.cptr(),v.ptr());
    } 
    else MultMV<false>(T(1),m.Transpose(),v,v); 
    return v; 
  }

  template <class T, int N, StorageType S, IndexStyle I2> 
  inline const VectorView<CT>& operator*=(
      const VectorView<CT>& v, const SmallMatrix<T,N,N,S,I2>& m)
  {
    TMVAssert(v.size() == N);
    if (v.step() == 1) { // conj is ok
      SmallVector<CT,N> v_copy(v);
      MultVM_1<N,N,S>(v_copy.cptr(),m.cptr(),v.ptr());
    } 
    else MultMV<false>(T(1),m.Transpose(),v,v); 
    return v; 
  }

  template <class T, int N, IndexStyle I>
  inline SmallVector<T,N,I>& operator*=(
      SmallVector<T,N,I>& v, const GenMatrix<T>& m)
  {
    TMVAssert(m.colsize() == N);
    TMVAssert(m.rowsize() == N);
    MultMV<false>(T(1),m.Transpose(),v,v); 
    return v; 
  }

  template <class T, int N, IndexStyle I>
  inline SmallVector<CT,N,I>& operator*=(
      SmallVector<CT,N,I>& v, const GenMatrix<T>& m)
  {
    TMVAssert(m.colsize() == N);
    TMVAssert(m.rowsize() == N);
    MultMV<false>(T(1),m.Transpose(),v,v); 
    return v; 
  }

  // v *= xm
  template <class T, class T1, int N, StorageType S, IndexStyle I1, IndexStyle I2>
  inline SmallVector<T,N,I1>& operator*=(
      SmallVector<T,N,I1>& v, const ProdXm<T,T1,N,N,S,I2>& pxm)
  {
    SmallVector<T,N> v_copy(v);
    MultVM<N,N,S>(pxm.GetX(),v_copy.cptr(),pxm.GetM().cptr(),v.ptr());
    return v; 
  }

  template <class T, int N, StorageType S, IndexStyle I1, IndexStyle I2> 
  inline SmallVector<CT,N,I1>& operator*=(
      SmallVector<CT,N,I1>& v, const ProdXm<T,T,N,N,S,I2>& pxm)
  {
    SmallVector<CT,N> v_copy(v);
    MultVM<N,N,S>(pxm.GetX(),v_copy.cptr(),pxm.GetM().cptr(),v.ptr());
    return v; 
  }

  template <class T, class T1, int N, StorageType S, IndexStyle I2>
  inline const VectorView<T>& operator*=(
      const VectorView<T>& v, const ProdXm<T,T1,N,N,S,I2>& pxm)
  {
    TMVAssert(v.size() == N);
    if (v.step() == 1 && !v.isconj()) {
      SmallVector<T,N> v_copy(v);
      MultVM<N,N,S>(pxm.GetX(),v_copy.cptr(),pxm.GetM().cptr(),v.ptr());
    } 
    else MultMV<false>(pxm.GetX(),pxm.GetM().Transpose(),v,v);
    return v;
  }

  template <class T, int N, StorageType S, IndexStyle I2> 
  inline const VectorView<CT>& operator*=(
      const VectorView<CT>& v, const ProdXm<T,T,N,N,S,I2>& pxm)
  {
    TMVAssert(v.size() == N);
    if (v.step() == 1) { // conj is ok
      SmallVector<CT,N> v_copy(v);
      MultVM<N,N,S>(pxm.GetX(),v_copy.cptr(),pxm.GetM().cptr(),v.ptr());
    } 
    else MultMV<false>(pxm.GetX(),pxm.GetM().Transpose(),v,v);
    return v;
  }

  template <class T, class T1, int N, StorageType S, IndexStyle I1>
  inline SmallVector<T,N,I1>& operator*=(
      SmallVector<T,N,I1>& v, const ProdXM<T,T1>& pxm)
  { 
    TMVAssert(pxm.colsize() == N);
    TMVAssert(pxm.rowsize() == N);
    MultMV<false>(pxm.GetX(),pxm.GetM().Transpose(),v,v);
    return v;
  }

  template <class T, int N, StorageType S, IndexStyle I1>
  inline SmallVector<CT,N,I1>& operator*=(
      SmallVector<CT,N,I1>& v, const ProdXM<T,T>& pxm)
  {
    TMVAssert(pxm.colsize() == N);
    TMVAssert(pxm.rowsize() == N);
    MultMV<false>(pxm.GetX(),pxm.GetM().Transpose(),v,v);
    return v;
  }

  // v += mv
  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I1, IndexStyle I2, IndexStyle I3>
  inline SmallVector<T,M,I1>& operator+=(
      SmallVector<T,M,I1>& v, const Prodmv_1<T,T1,T2,M,N,S,I2,I3>& pmv)
  {
    AddMultMV_1<M,N,S>(pmv.GetM().cptr(),pmv.GetV().cptr(),v.ptr());
    return v;
  }

  template <class T, int M, int N, StorageType S, IndexStyle I1, IndexStyle I2, IndexStyle I3>
  inline SmallVector<CT,M,I1>& operator+=(
      SmallVector<CT,M,I1>& v, const Prodmv_1<T,T,T,M,N,S,I2,I3>& pmv)
  {
    AddMultMV_1<M,N,S>(pmv.GetM().cptr(),pmv.GetV().cptr(),v.ptr());
    return v;
  }

  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I2, IndexStyle I3>
  inline const VectorView<T>& operator+=(
      const VectorView<T>& v, const Prodmv_1<T,T1,T2,M,N,S,I2,I3>& pmv)
  {
    TMVAssert(v.size() == M);
    if (v.step() == 1 && !v.isconj()) 
      AddMultMV_1<M,N,S>(pmv.GetM().cptr(),pmv.GetV().cptr(),v.ptr());
    else 
      MultMV<true>(T(1),pmv.GetM().View(),pmv.GetV().View(),v); 
    return v; 
  }

  template <class T, int M, int N, StorageType S, IndexStyle I2, IndexStyle I3>
  inline const VectorView<CT>& operator+=(
      const VectorView<CT>& v, const Prodmv_1<T,T,T,M,N,S,I2,I3>& pmv)
  {
    TMVAssert(v.size() == M);
    if (v.step() == 1)  // conj is ok
      AddMultMV_1<M,N,S>(pmv.GetM().cptr(),pmv.GetV().cptr(),v.ptr());
    else 
      MultMV<true>(T(1),pmv.GetM().View(),pmv.GetV().View(),v); 
    return v; 
  }

  // v -= mv
  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I1, IndexStyle I2, IndexStyle I3>
  inline SmallVector<T,M,I1>& operator-=(
      SmallVector<T,M,I1>& v, const Prodmv_1<T,T1,T2,M,N,S,I2,I3>& pmv)
  { 
    AddMultMV_m1<M,N,S>(pmv.GetM().cptr(), pmv.GetV().cptr(), v.ptr());
    return v;
  }

  template <class T, int M, int N, StorageType S, IndexStyle I1, IndexStyle I2, IndexStyle I3>
  inline SmallVector<CT,M,I1>& operator-=(
      SmallVector<CT,M,I1>& v, const Prodmv_1<T,T,T,M,N,S,I2,I3>& pmv)
  {
    AddMultMV_m1<M,N,S>(pmv.GetM().cptr(),pmv.GetV().cptr(),v.ptr());
    return v;
  }

  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I2, IndexStyle I3>
  inline const VectorView<T>& operator-=(
      const VectorView<T>& v, const Prodmv_1<T,T1,T2,M,N,S,I2,I3>& pmv)
  {
    TMVAssert(v.size() == M);
    if (v.step() == 1 && !v.isconj()) 
      AddMultMV_m1<M,N,S>(pmv.GetM().cptr(),pmv.GetV().cptr(),v.ptr());
    else 
      MultMV<true>(T(-1),pmv.GetM().View(),pmv.GetV().View(),v); 
    return v; 
  }

  template <class T, int M, int N, StorageType S, IndexStyle I2, IndexStyle I3>
  inline const VectorView<CT>& operator-=(
      const VectorView<CT>& v, const Prodmv_1<T,T,T,M,N,S,I2,I3>& pmv)
  {
    TMVAssert(v.size() == M);
    if (v.step() == 1)  // conj is ok
      AddMultMV_m1<M,N,S>(pmv.GetM().cptr(),pmv.GetV().cptr(),v.ptr());
    else 
      MultMV<true>(T(-1),pmv.GetM().View(),pmv.GetV().View(),v); 
    return v; 
  }

  // v += vm
  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I1, IndexStyle I2, IndexStyle I3>
  inline SmallVector<T,N,I1>& operator+=(
      SmallVector<T,N,I1>& v, const Prodvm_1<T,T1,T2,M,N,S,I2,I3>& pvm)
  { 
    AddMultVM_1<M,N,S>(pvm.GetV().cptr(),pvm.GetM().cptr(),v.ptr());
    return v;
  }

  template <class T, int M, int N, StorageType S, IndexStyle I1, IndexStyle I2, IndexStyle I3>
  inline SmallVector<CT,N,I1>& operator+=(
      SmallVector<CT,N,I1>& v, const Prodvm_1<T,T,T,M,N,S,I2,I3>& pvm)
  {
    AddMultVM_1<M,N,S>(pvm.GetV().cptr(),pvm.GetM().cptr(),v.ptr());
    return v;
  }

  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I2, IndexStyle I3>
  inline const VectorView<T>& operator+=(
      const VectorView<T>& v, const Prodvm_1<T,T1,T2,M,N,S,I2,I3>& pvm)
  {
    TMVAssert(v.size() == N);
    if (v.step() == 1 && !v.isconj()) 
      AddMultVM_1<M,N,S>(pvm.GetV().cptr(),pvm.GetM().cptr(),v.ptr());
    else 
      MultMV<true>(T(1),pvm.GetM().Transpose(),pvm.GetV().View(),v); 
    return v; 
  }

  template <class T, int M, int N, StorageType S, IndexStyle I2, IndexStyle I3>
  inline const VectorView<CT>& operator+=(
      const VectorView<CT>& v, const Prodvm_1<T,T,T,M,N,S,I2,I3>& pvm)
  {
    TMVAssert(v.size() == N);
    if (v.step() == 1)  // conj is ok
      AddMultVM_1<M,N,S>(pvm.GetV().cptr(),pvm.GetM().cptr(),v.ptr());
    else 
      MultMV<true>(T(1),pvm.GetM().Transpose(),pvm.GetV().View(),v); 
    return v; 
  }

  // v -= vm
  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I1, IndexStyle I2, IndexStyle I3>
  inline SmallVector<T,N,I1>& operator-=(
      SmallVector<T,N,I1>& v, const Prodvm_1<T,T1,T2,M,N,S,I2,I3>& pvm)
  {
    AddMultVM_m1<M,N,S>(pvm.GetV().cptr(),pvm.GetM().cptr(),v.ptr());
    return v;
  }

  template <class T, int M, int N, StorageType S, IndexStyle I1, IndexStyle I2, IndexStyle I3>
  inline SmallVector<CT,N,I1>& operator-=(
      SmallVector<CT,N,I1>& v, const Prodvm_1<T,T,T,M,N,S,I2,I3>& pvm)
  {
    AddMultVM_m1<M,N,S>(pvm.GetV().cptr(),pvm.GetM().cptr(),v.ptr());
    return v;
  }

  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I2, IndexStyle I3>
  inline const VectorView<T>& operator-=(
      const VectorView<T>& v, const Prodvm_1<T,T1,T2,M,N,S,I2,I3>& pvm)
  {
    TMVAssert(v.size() == N);
    if (v.step() == 1 && !v.isconj()) 
      AddMultVM_m1<M,N,S>(pvm.GetV().cptr(),pvm.GetM().cptr(),v.ptr());
    else 
      MultMV<true>(T(-1),pvm.GetM().Transpose(),pvm.GetV().View(),v); 
    return v; 
  }

  template <class T, int M, int N, StorageType S, IndexStyle I2, IndexStyle I3>
  inline const VectorView<CT>& operator-=(
      const VectorView<CT>& v, const Prodvm_1<T,T,T,M,N,S,I2,I3>& pvm)
  {
    TMVAssert(v.size() == N);
    if (v.step() == 1)  // conj is ok
      AddMultVM_m1<M,N,S>(pvm.GetV().cptr(),pvm.GetM().cptr(),v.ptr());
    else 
      MultMV<true>(T(-1),pvm.GetM().Transpose(),pvm.GetV().View(),v); 
    return v; 
  }

  // v += xmv
  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I1, IndexStyle I2, IndexStyle I3>
  inline SmallVector<T,M,I1>& operator+=(
      SmallVector<T,M,I1>& v, const Prodmv<T,T1,T2,M,N,S,I2,I3>& pmv)
  {
    AddMultMV<M,N,S>(pmv.GetX(),pmv.GetM().cptr(),pmv.GetV().cptr(),v.ptr());
    return v;
  }

  template <class T, int M, int N, StorageType S, IndexStyle I1, IndexStyle I2, IndexStyle I3>
  inline SmallVector<CT,M,I1>& operator+=(
      SmallVector<CT,M,I1>& v, const Prodmv<T,T,T,M,N,S,I2,I3>& pmv)
  {
    AddMultMV<M,N,S>(pmv.GetX(),pmv.GetM().cptr(),pmv.GetV().cptr(),v.ptr());
    return v;
  }

  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I2, IndexStyle I3>
  inline const VectorView<T>& operator+=(
      const VectorView<T>& v, const Prodmv<T,T1,T2,M,N,S,I2,I3>& pmv)
  {
    TMVAssert(v.size() == M);
    if (v.step() == 1 && !v.isconj()) 
      AddMultMV<M,N,S>(pmv.GetX(),pmv.GetM().cptr(),pmv.GetV().cptr(),
          v.ptr());
    else 
      MultMV<true>(pmv.GetX(),pmv.GetM().View(),pmv.GetV().View(),v);
    return v; 
  }

  template <class T, int M, int N, StorageType S, IndexStyle I2, IndexStyle I3>
  inline const VectorView<CT>& operator+=(
      const VectorView<CT>& v, const Prodmv<T,T,T,M,N,S,I2,I3>& pmv)
  {
    TMVAssert(v.size() == M);
    if (v.step() == 1)  // conj is ok
      AddMultMV<M,N,S>(pmv.GetX(),pmv.GetM().cptr(),pmv.GetV().cptr(),
          v.ptr());
    else 
      MultMV<true>(pmv.GetX(),pmv.GetM().View(),pmv.GetV().View(),v);
    return v; 
  }

  template <class T, class T1, class T2, int M, IndexStyle I1>
  inline SmallVector<T,M,I1>& operator+=(
      SmallVector<T,M,I1>& v, const ProdMV<T,T1,T2>& pmv)
  {
    TMVAssert(pmv.size() == M);
    MultMV<true>(pmv.GetX(),pmv.GetM().View(),pmv.GetV().View(),v.View()); 
    return v; 
  }

  template <class T, int M, IndexStyle I1> 
  inline SmallVector<CT,M,I1>& operator+=(
      SmallVector<CT,M,I1>& v, const ProdMV<T,T,T>& pmv)
  {
    TMVAssert(pmv.size() == M);
    MultMV<true>(pmv.GetX(),pmv.GetM().View(),pmv.GetV().View(),v.View()); 
    return v; 
  }

  // v -= xmv
  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I1, IndexStyle I2, IndexStyle I3>
  inline SmallVector<T,M,I1>& operator-=(
      SmallVector<T,M,I1>& v, const Prodmv<T,T1,T2,M,N,S,I2,I3>& pmv)
  {
    AddMultMV<M,N,S>(-pmv.GetX(),pmv.GetM().cptr(),pmv.GetV().cptr(),v.ptr());
    return v;
  }

  template <class T, int M, int N, StorageType S, IndexStyle I1, IndexStyle I2, IndexStyle I3>
  inline SmallVector<CT,M,I1>& operator-=(
      SmallVector<CT,M,I1>& v, const Prodmv<T,T,T,M,N,S,I2,I3>& pmv)
  {
    AddMultMV<M,N,S>(-pmv.GetX(),pmv.GetM().cptr(),pmv.GetV().cptr(),v.ptr());
    return v;
  }

  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I2, IndexStyle I3>
  inline const VectorView<T>& operator-=(
      const VectorView<T>& v, const Prodmv<T,T1,T2,M,N,S,I2,I3>& pmv)
  {
    TMVAssert(v.size() == M);
    if (v.step() == 1 && !v.isconj()) 
      AddMultMV<M,N,S>(-pmv.GetX(),pmv.GetM().cptr(),pmv.GetV().cptr(),
          v.ptr());
    else 
      MultMV<true>(-pmv.GetX(),pmv.GetM().View(),pmv.GetV().View(),v);
    return v; 
  }

  template <class T, int M, int N, StorageType S, IndexStyle I2, IndexStyle I3>
  inline const VectorView<CT>& operator-=(
      const VectorView<CT>& v, const Prodmv<T,T,T,M,N,S,I2,I3>& pmv)
  {
    TMVAssert(v.size() == M);
    if (v.step() == 1)  // conj is ok
      AddMultMV<M,N,S>(-pmv.GetX(),pmv.GetM().cptr(),pmv.GetV().cptr(),
          v.ptr());
    else 
      MultMV<true>(-pmv.GetX(),pmv.GetM().View(),pmv.GetV().View(),v);
    return v; 
  }

  template <class T, class T1, class T2, int M, IndexStyle I1>
  inline SmallVector<T,M,I1>& operator-=(
      SmallVector<T,M,I1>& v, const ProdMV<T,T1,T2>& pmv)
  {
    TMVAssert(pmv.size() == M);
    MultMV<true>(-pmv.GetX(),pmv.GetM().View(),pmv.GetV().View(),v.View()); 
    return v; 
  }

  template <class T, int M, IndexStyle I1> 
  inline SmallVector<CT,M,I1>& operator-=(
      SmallVector<CT,M,I1>& v, const ProdMV<T,T,T>& pmv)
  {
    TMVAssert(pmv.size() == M);
    MultMV<true>(-pmv.GetX(),pmv.GetM().View(),pmv.GetV().View(),v.View()); 
    return v; 
  }

  // v += xvm
  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I1, IndexStyle I2, IndexStyle I3>
  inline SmallVector<T,N,I1>& operator+=(
      SmallVector<T,N,I1>& v, const Prodvm<T,T1,T2,M,N,S,I2,I3>& pvm)
  {
    AddMultVM<M,N,S>(pvm.GetX(),pvm.GetV().cptr(),pvm.GetM().cptr(),v.ptr());
    return v;
  }

  template <class T, int M, int N, StorageType S, IndexStyle I1, IndexStyle I2, IndexStyle I3>
  inline SmallVector<CT,N,I1>& operator+=(
      SmallVector<CT,N,I1>& v, const Prodvm<T,T,T,M,N,S,I2,I3>& pvm)
  {
    AddMultVM<M,N,S>(pvm.GetX(),pvm.GetV().cptr(),pvm.GetM().cptr(),v.ptr());
    return v;
  }

  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I2, IndexStyle I3>
  inline const VectorView<T>& operator+=(
      const VectorView<T>& v, const Prodvm<T,T1,T2,M,N,S,I2,I3>& pvm)
  {
    TMVAssert(v.size() == N);
    if (v.step() == 1 && !v.isconj()) 
      AddMultVM<M,N,S>(pvm.GetX(),pvm.GetV().cptr(),pvm.GetM().cptr(),
          v.ptr());
    else 
      MultMV<true>(pvm.GetX(),pvm.GetM().Transpose(),pvm.GetV().View(),v);
    return v; 
  }

  template <class T, int M, int N, StorageType S, IndexStyle I2, IndexStyle I3>
  inline const VectorView<CT>& operator+=(
      const VectorView<CT>& v, const Prodvm<T,T,T,M,N,S,I2,I3>& pvm)
  {
    TMVAssert(v.size() == N);
    if (v.step() == 1)  // conj is ok
      AddMultVM<M,N,S>(pvm.GetX(),pvm.GetV().cptr(),pvm.GetM().cptr(),
          v.ptr());
    else 
      MultMV<true>(pvm.GetX(),pvm.GetM().Transpose(),pvm.GetV().View(),v);
    return v; 
  }

  // v -= xvm
  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I1, IndexStyle I2, IndexStyle I3>
  inline SmallVector<T,N,I1>& operator-=(
      SmallVector<T,N,I1>& v, const Prodvm<T,T1,T2,M,N,S,I2,I3>& pvm)
  {
    AddMultVM<M,N,S>(-pvm.GetX(),pvm.GetV().cptr(),pvm.GetM().cptr(),v.ptr());
    return v;
  }

  template <class T, int M, int N, StorageType S, IndexStyle I1, IndexStyle I2, IndexStyle I3>
  inline SmallVector<CT,N,I1>& operator-=(
      SmallVector<CT,N,I1>& v, const Prodvm<T,T,T,M,N,S,I2,I3>& pvm)
  {
    AddMultVM<M,N,S>(-pvm.GetX(),pvm.GetV().cptr(),pvm.GetM().cptr(),v.ptr());
    return v;
  }

  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I2, IndexStyle I3>
  inline const VectorView<T>& operator-=(
      const VectorView<T>& v, const Prodvm<T,T1,T2,M,N,S,I2,I3>& pvm)
  {
    TMVAssert(v.size() == N);
    if (v.step() == 1 && !v.isconj()) 
      AddMultVM<M,N,S>(-pvm.GetX(),pvm.GetV().cptr(),pvm.GetM().cptr(),
          v.ptr());
    else 
      MultMV<true>(-pvm.GetX(),pvm.GetM().Transpose(),pvm.GetV().View(),v);
    return v; 
  }

  template <class T, int M, int N, StorageType S, IndexStyle I2, IndexStyle I3>
  inline const VectorView<CT>& operator-=(
      const VectorView<CT>& v, const Prodvm<T,T,T,M,N,S,I2,I3>& pvm)
  {
    TMVAssert(v.size() == N);
    if (v.step() == 1)  // conj is ok
      AddMultVM<M,N,S>(-pvm.GetX(),pvm.GetV().cptr(),pvm.GetM().cptr(),
          v.ptr());
    else 
      MultMV<true>(-pvm.GetX(),pvm.GetM().Transpose(),pvm.GetV().View(),v);
    return v; 
  }

#define GENMATRIX1 SmallMatrix
#define GENMATRIX2 SmallVector
#define PRODMM Prodmv
#define PRODXM1 ProdXm
#define PRODXM2 ProdXv
#define PRODMM_1 Prodmv_1
#define X1 ,M,N,S,I1
#define X2 ,N,I2
#define X3 ,M,N,S,I1,I2
#define Y ,int M,int N,StorageType S, IndexStyle I1, IndexStyle I2
#define GETM1 .GetM()
#define GETM2 .GetV()
#include "tmv/TMV_AuxProdMM.h"
#define PRODMM_1 Prodmv_1
#define X3 ,M,N,S,I1,I2
#define Y ,int M,int N,StorageType S, IndexStyle I1, IndexStyle I2
#define GETM1 .GetM()
#define GETM2 .GetV()
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 SmallVector
#define GENMATRIX2 SmallMatrix
#define PRODMM Prodvm
#define PRODXM1 ProdXv
#define PRODXM2 ProdXm
#define PRODMM_1 Prodvm_1
#define X1 ,M,I1
#define X2 ,M,N,S,I2
#define X3 ,M,N,S,I1,I2
#define Y ,int M,int N,StorageType S, IndexStyle I1, IndexStyle I2
#define GETM1 .GetV()
#define GETM2 .GetM()
#include "tmv/TMV_AuxProdMM.h"
#define PRODMM_1 Prodvm_1
#define X3 ,M,N,S,I1,I2
#define Y ,int M,int N,StorageType S, IndexStyle I1, IndexStyle I2
#define GETM1 .GetV()
#define GETM2 .GetM()
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 SmallMatrix
#define GENMATRIX2 GenVector
#define PRODMM ProdmV
#define PRODXM1 ProdXm
#define PRODXM2 ProdXV
#define X1 ,M,N,S,I1
#define X3 ,M,N,S,I1
#define Y ,int M,int N,StorageType S, IndexStyle I1
#define GETM1 .GetM()
#define GETM2 .GetV()
#include "tmv/TMV_AuxProdMM.h"
#define X3 ,M,N,S,I1
#define Y ,int M,int N,StorageType S, IndexStyle I1
#define GETM1 .GetM()
#define GETM2 .GetV()
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 SmallVector
#define PRODMM ProdMv
#define PRODXM1 ProdXM
#define PRODXM2 ProdXv
#define X2 ,N,I2
#define X3 ,N,I2
#define Y ,int N, IndexStyle I2
#define GETM1 .GetM()
#define GETM2 .GetV()
#include "tmv/TMV_AuxProdMM.h"
#define X3 ,N,I2
#define Y ,int N, IndexStyle I2
#define GETM1 .GetM()
#define GETM2 .GetV()
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 SmallVector
#define GENMATRIX2 GenMatrix
#define PRODMM ProdvM
#define PRODXM1 ProdXv
#define PRODXM2 ProdXM
#define X1 ,M,I1
#define X3 ,M,I1
#define Y ,int M, IndexStyle I1
#define GETM1 .GetV()
#define GETM2 .GetM()
#include "tmv/TMV_AuxProdMM.h"
#define X3 ,M,I1
#define Y ,int M, IndexStyle I1
#define GETM1 .GetV()
#define GETM2 .GetM()
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX1 GenVector
#define GENMATRIX2 SmallMatrix
#define PRODMM ProdVm
#define PRODXM1 ProdXV
#define PRODXM2 ProdXm
#define X2 ,M,N,S,I2
#define X3 ,M,N,S,I2
#define Y ,int M,int N,StorageType S, IndexStyle I2
#define GETM1 .GetV()
#define GETM2 .GetM()
#include "tmv/TMV_AuxProdMM.h"
#define X3 ,M,N,S,I2
#define Y ,int M,int N,StorageType S, IndexStyle I2
#define GETM1 .GetV()
#define GETM2 .GetM()
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODMM
#undef PRODXM1
#undef PRODXM2


  //
  // Scalar / Matrix
  //

  template <class T, class Tm, int M, int N, StorageType S, IndexStyle I> 
  class QuotXm_1 : 
    public SmallMatrixComposite<T,N,M> 
  {
  public:
    inline QuotXm_1(const T DEBUGPARAM(_x), 
        const SmallMatrix<Tm,M,N,S,I>& _m) :
      m(_m) { TMVAssert(_x == T(1)); }
    inline const SmallMatrix<Tm,M,N,S,I>& GetM() const { return m; }
    inline void AssignTom(
        SmallMatrix<RealType(T),N,M,RowMajor,CStyle>& m0) const
    { TMVAssert(IsReal(T())); DoInverse(m,m0); }
    inline void AssignTom(
        SmallMatrix<RealType(T),N,M,ColMajor,CStyle>& m0) const
    { TMVAssert(IsReal(T())); DoInverse(m,m0); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),N,M,RowMajor,CStyle>& m0) const
    { DoInverse(m,m0); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),N,M,ColMajor,CStyle>& m0) const
    { DoInverse(m,m0); }
    inline void AssignTom(
        SmallMatrix<RealType(T),N,M,RowMajor,FortranStyle>& m0) const
    { TMVAssert(IsReal(T())); DoInverse(m,m0); }
    inline void AssignTom(
        SmallMatrix<RealType(T),N,M,ColMajor,FortranStyle>& m0) const
    { TMVAssert(IsReal(T())); DoInverse(m,m0); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),N,M,RowMajor,FortranStyle>& m0) const
    { DoInverse(m,m0); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),N,M,ColMajor,FortranStyle>& m0) const
    { DoInverse(m,m0); }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    { TMVAssert(IsReal(T())); m.View().Inverse(m0); }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    { m.View().Inverse(m0); }
  private:
    const SmallMatrix<Tm,M,N,S,I>& m;
  };

  template <class T, class Tm, int M, int N, StorageType S, IndexStyle I> 
  class QuotXm : 
    public SmallMatrixComposite<T,N,M> 
  {
  public:
    inline QuotXm(const T _x, const SmallMatrix<Tm,M,N,S,I>& _m) : 
      x(_x), m(_m) {}
    inline StorageType stor() const { return S; }
    inline T GetX() const { return x; }
    inline const SmallMatrix<Tm,M,N,S,I>& GetM() const { return m; }
    inline void AssignTom(
        SmallMatrix<RealType(T),N,M,RowMajor,CStyle>& m0) const
    { TMVAssert(IsReal(T())); DoInverse(m,m0); MultXV<M*N>(x,m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<RealType(T),N,M,ColMajor,CStyle>& m0) const
    { TMVAssert(IsReal(T())); DoInverse(m,m0); MultXV<M*N>(x,m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),N,M,RowMajor,CStyle>& m0) const
    { DoInverse(m,m0); MultXV<M*N>(x,m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),N,M,ColMajor,CStyle>& m0) const
    { DoInverse(m,m0); MultXV<M*N>(x,m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<RealType(T),N,M,RowMajor,FortranStyle>& m0) const
    { TMVAssert(IsReal(T())); DoInverse(m,m0); MultXV<M*N>(x,m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<RealType(T),N,M,ColMajor,FortranStyle>& m0) const
    { TMVAssert(IsReal(T())); DoInverse(m,m0); MultXV<M*N>(x,m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),N,M,RowMajor,FortranStyle>& m0) const
    { DoInverse(m,m0); MultXV<M*N>(x,m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),N,M,ColMajor,FortranStyle>& m0) const
    { DoInverse(m,m0); MultXV<M*N>(x,m0.ptr()); }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    { TMVAssert(IsReal(T())); m.View().Inverse(m0); MultXM(x,m0); }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    { m.View().Inverse(m0); MultXM(x,m0); }
  private:
    const T x;
    const SmallMatrix<Tm,M,N,S,I>& m;
  };

#define GENMATRIX SmallMatrix
#define PRODXM ProdXm
#define QUOTXM QuotXm
#define X ,M,N,S,I
#define Y ,int M,int N,StorageType S, IndexStyle I
#include "tmv/TMV_AuxQuotXM.h"
#define X ,M,N,S,I
#define Y ,int M,int N,StorageType S, IndexStyle I
#define QUOTXM_1 QuotXm_1
#include "tmv/TMV_AuxQuotXMa.h"
#undef GENMATRIX
#undef PRODXM
#undef QUOTXM


  //
  // Vector / % Matrix 
  //

  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I1, IndexStyle I2> 
  class Quotvm_1 :
    public SmallVectorComposite<T,N>
  {
  public:
    inline Quotvm_1(
        const T DEBUGPARAM(_x), const SmallVector<T1,M,I1>& _v,
        const SmallMatrix<T2,M,N,S,I2>& _m) :
      v(_v), m(_m) { TMVAssert(_x == T(1)); }
    inline const SmallVector<T1,M,I1>& GetV() const { return v; }
    inline const SmallMatrix<T2,M,N,S,I2>& GetM() const { return m; }
    inline void AssignTov(SmallVector<RealType(T),N,CStyle>& v0) const
    { TMVAssert(IsReal(T())); DoLDiv(m,v,v0); }
    inline void AssignTov(SmallVector<ComplexType(T),N,CStyle>& v0) const
    { DoLDiv(m,v,v0); }
    inline void AssignTov(SmallVector<RealType(T),N,FortranStyle>& v0) const
    { TMVAssert(IsReal(T())); DoLDiv(m,v,v0); }
    inline void AssignTov(SmallVector<ComplexType(T),N,FortranStyle>& v0) const
    { DoLDiv(m,v,v0); }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    { 
      TMVAssert(v0.size() == N);
      TMVAssert(IsReal(T())); 
      m.View().LDiv(v.View(),v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    { 
      TMVAssert(v0.size() == N);
      m.View().LDiv(v.View(),v0);
    }
  private:
    const SmallVector<T1,M,I1>& v;
    const SmallMatrix<T2,M,N,S,I2>& m;
  };

  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I1, IndexStyle I2> 
  class Quotvm :
    public SmallVectorComposite<T,N>
  {
  public:
    inline Quotvm(const T _x, const SmallVector<T1,M,I1>& _v,
        const SmallMatrix<T2,M,N,S,I2>& _m) :
      x(_x), v(_v), m(_m) {}
    inline T GetX() const { return x; }
    inline const SmallVector<T1,M,I1>& GetV() const { return v; }
    inline const SmallMatrix<T2,M,N,S,I2>& GetM() const { return m; }
    inline void AssignTov(SmallVector<RealType(T),N,CStyle>& v0) const
    { TMVAssert(IsReal(T())); DoLDiv(m,v,v0); MultXV<N>(x,v0.ptr()); }
    inline void AssignTov(SmallVector<ComplexType(T),N,CStyle>& v0) const
    { DoLDiv(m,v,v0); MultXV<N>(x,v0.ptr()); }
    inline void AssignTov(SmallVector<RealType(T),N,FortranStyle>& v0) const
    { TMVAssert(IsReal(T())); DoLDiv(m,v,v0); MultXV<N>(x,v0.ptr()); }
    inline void AssignTov(SmallVector<ComplexType(T),N,FortranStyle>& v0) const
    { DoLDiv(m,v,v0); MultXV<N>(x,v0.ptr()); }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    { 
      TMVAssert(v0.size() == N);
      TMVAssert(IsReal(T())); 
      m.View().LDiv(v.View(),v0);
      MultXV(x,v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    { 
      TMVAssert(v0.size() == N);
      m.View().LDiv(v.View(),v0);
      MultXV(x,v0);
    }
  private:
    const T x;
    const SmallVector<T1,M,I1>& v;
    const SmallMatrix<T2,M,N,S,I2>& m;
  };

  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I> 
  class QuotVm_1 :
    public VectorComposite<T>
  {
  public:
    inline QuotVm_1(const T DEBUGPARAM(_x), const GenVector<T1>& _v,
        const SmallMatrix<T2,M,N,S,I>& _m) :
      v(_v), m(_m) { TMVAssert(_x == T(1)); TMVAssert(v.size() == M); }
    inline size_t size() const { return N; }
    inline const GenVector<T1>& GetV() const { return v; }
    inline const SmallMatrix<T2,M,N,S,I>& GetM() const { return m; }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    { 
      TMVAssert(v0.size() == size());
      TMVAssert(IsReal(T())); 
      m.View().LDiv(v,v0);
      MultXV(T(1),v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    { 
      TMVAssert(v0.size() == size());
      m.View().LDiv(v,v0);
      MultXV(T(1),v0);
    }
  private:
    const GenVector<T1>& v;
    const SmallMatrix<T2,M,N,S,I>& m;
  };

  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I> 
  class QuotVm :
    public VectorComposite<T>
  {
  public:
    inline QuotVm(const T _x, const GenVector<T1>& _v,
        const SmallMatrix<T2,M,N,S,I>& _m) :
      x(_x), v(_v), m(_m) { TMVAssert(v.size() == M); }
    inline size_t size() const { return N; }
    inline T GetX() const { return x; }
    inline const GenVector<T1>& GetV() const { return v; }
    inline const SmallMatrix<T2,M,N,S,I>& GetM() const { return m; }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    { 
      TMVAssert(v0.size() == size());
      TMVAssert(IsReal(T())); 
      m.View().LDiv(v,v0);
      MultXV(x,v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    { 
      TMVAssert(v0.size() == size());
      m.View().LDiv(v,v0);
      MultXV(x,v0);
    }
  private:
    const T x;
    const GenVector<T1>& v;
    const SmallMatrix<T2,M,N,S,I>& m;
  };

  template <class T, class T1, class T2, int M, IndexStyle I> 
  class QuotvM :
    public VectorComposite<T>
  {
  public:
    inline QuotvM(const T _x, const SmallVector<T1,M,I>& _v,
        const GenMatrix<T2>& _m) :
      x(_x), v(_v), m(_m) { TMVAssert(m.colsize() == M); }
    inline size_t size() const { return m.rowsize(); }
    inline T GetX() const { return x; }
    inline const SmallVector<T1,M,I>& GetV() const { return v; }
    inline const GenMatrix<T2>& GetM() const { return m; }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    { 
      TMVAssert(v0.size() == size());
      TMVAssert(IsReal(T())); 
      m.LDiv(v.View(),v0);
      MultXV(x,v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    { 
      TMVAssert(v0.size() == size());
      m.LDiv(v.View(),v0);
      MultXV(x,v0);
    }
  private:
    const T x;
    const SmallVector<T1,M,I>& v;
    const GenMatrix<T2>& m;
  };

  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I1, IndexStyle I2> 
  class RQuotvm_1 :
    public SmallVectorComposite<T,M>
  {
  public:
    inline RQuotvm_1(
        const T DEBUGPARAM(_x), const SmallVector<T1,N,I1>& _v,
        const SmallMatrix<T2,M,N,S,I2>& _m) :
      v(_v), m(_m) { TMVAssert(_x == T(1)); }
    inline const SmallVector<T1,N,I1>& GetV() const { return v; }
    inline const SmallMatrix<T2,M,N,S,I2>& GetM() const { return m; }
    inline void AssignTov(SmallVector<RealType(T),M,CStyle>& v0) const
    { TMVAssert(IsReal(T())); DoRDiv(m,v,v0); }
    inline void AssignTov(SmallVector<ComplexType(T),M,CStyle>& v0) const
    { DoRDiv(m,v,v0); }
    inline void AssignTov(SmallVector<RealType(T),M,FortranStyle>& v0) const
    { TMVAssert(IsReal(T())); DoRDiv(m,v,v0); }
    inline void AssignTov(SmallVector<ComplexType(T),M,FortranStyle>& v0) const
    { DoRDiv(m,v,v0); }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    { 
      TMVAssert(v0.size() == M);
      TMVAssert(IsReal(T())); 
      m.View().RDiv(v.View(),v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    { 
      TMVAssert(v0.size() == M);
      m.View().RDiv(v.View(),v0);
    }
  private:
    const SmallVector<T1,N,I1>& v;
    const SmallMatrix<T2,M,N,S,I2>& m;
  };

  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I1, IndexStyle I2> 
  class RQuotvm :
    public SmallVectorComposite<T,M>
  {
  public:
    inline RQuotvm(const T _x, const SmallVector<T1,N,I1>& _v,
        const SmallMatrix<T2,M,N,S,I2>& _m) :
      x(_x), v(_v), m(_m) {}
    inline T GetX() const { return x; }
    inline const SmallVector<T1,N,I1>& GetV() const { return v; }
    inline const SmallMatrix<T2,M,N,S,I2>& GetM() const { return m; }
    inline void AssignTov(SmallVector<RealType(T),M,CStyle>& v0) const
    { TMVAssert(IsReal(T())); DoRDiv(m,v,v0); MultXV<M>(x,v0.ptr()); }
    inline void AssignTov(SmallVector<ComplexType(T),M,CStyle>& v0) const
    { DoRDiv(m,v,v0); MultXV<M>(x,v0.ptr()); }
    inline void AssignTov(SmallVector<RealType(T),M,FortranStyle>& v0) const
    { TMVAssert(IsReal(T())); DoRDiv(m,v,v0); MultXV<M>(x,v0.ptr()); }
    inline void AssignTov(SmallVector<ComplexType(T),M,FortranStyle>& v0) const
    { DoRDiv(m,v,v0); MultXV<M>(x,v0.ptr()); }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    { 
      TMVAssert(v0.size() == M);
      TMVAssert(IsReal(T())); 
      m.View().RDiv(v.View(),v0);
      MultXV(x,v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    { 
      TMVAssert(v0.size() == M);
      m.View().RDiv(v.View(),v0);
      MultXV(x,v0);
    }
  private:
    const T x;
    const SmallVector<T1,N,I1>& v;
    const SmallMatrix<T2,M,N,S,I2>& m;
  };

  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I> 
  class RQuotVm_1 :
    public VectorComposite<T>
  {
  public:
    inline RQuotVm_1(const T DEBUGPARAM(_x), const GenVector<T1>& _v,
        const SmallMatrix<T2,M,N,S,I>& _m) :
      v(_v), m(_m) { TMVAssert(_x==T(1)); TMVAssert(v.size() == N); }
    inline size_t size() const { return M; }
    inline const GenVector<T1>& GetV() const { return v; }
    inline const SmallMatrix<T2,M,N,S,I>& GetM() const { return m; }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    { 
      TMVAssert(v0.size() == size());
      TMVAssert(IsReal(T())); 
      m.View().RDiv(v,v0);
      MultXV(T(1),v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    { 
      TMVAssert(v0.size() == size());
      m.View().RDiv(v,v0);
      MultXV(T(1),v0);
    }
  private:
    const GenVector<T1>& v;
    const SmallMatrix<T2,M,N,S,I>& m;
  };

  template <class T, class T1, class T2, int M, int N, StorageType S, IndexStyle I> 
  class RQuotVm :
    public VectorComposite<T>
  {
  public:
    inline RQuotVm(const T _x, const GenVector<T1>& _v,
        const SmallMatrix<T2,M,N,S,I>& _m) :
      x(_x), v(_v), m(_m) { TMVAssert(v.size() == N); }
    inline size_t size() const { return M; }
    inline T GetX() const { return x; }
    inline const GenVector<T1>& GetV() const { return v; }
    inline const SmallMatrix<T2,M,N,S,I>& GetM() const { return m; }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    { 
      TMVAssert(v0.size() == size());
      TMVAssert(IsReal(T())); 
      m.View().RDiv(v,v0);
      MultXV(x,v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    { 
      TMVAssert(v0.size() == size());
      m.View().RDiv(v,v0);
      MultXV(x,v0);
    }
  private:
    const T x;
    const GenVector<T1>& v;
    const SmallMatrix<T2,M,N,S,I>& m;
  };

  template <class T, class T1, class T2, int N, IndexStyle I> 
  class RQuotvM :
    public VectorComposite<T>
  {
  public:
    inline RQuotvM(const T _x, const SmallVector<T1,N,I>& _v,
        const GenMatrix<T2>& _m) :
      x(_x), v(_v), m(_m) { TMVAssert(m.rowsize() == N); }
    inline size_t size() const { return m.colsize(); }
    inline T GetX() const { return x; }
    inline const SmallVector<T1,N,I>& GetV() const { return v; }
    inline const GenMatrix<T2>& GetM() const { return m; }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    { 
      TMVAssert(v0.size() == size());
      TMVAssert(IsReal(T())); 
      m.RDiv(v.View(),v0);
      MultXV(x,v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    { 
      TMVAssert(v0.size() == size());
      m.RDiv(v.View(),v0);
      MultXV(x,v0);
    }
  private:
    const T x;
    const SmallVector<T1,N,I>& v;
    const GenMatrix<T2>& m;
  };

  template <class T, int N, StorageType S, IndexStyle I1, IndexStyle I2> 
  inline SmallVector<T,N,I1>& operator/=(
      SmallVector<T,N,I1>& v, const SmallMatrix<T,N,N,S,I2>& m)
  { 
    DoLDivEq(m,v);
    return v; 
  }

  template <class T, int N, StorageType S, IndexStyle I1, IndexStyle I2> 
  inline SmallVector<CT,N,I1>& operator/=(
      SmallVector<CT,N,I1>& v, const SmallMatrix<T,N,N,S,I2>& m)
  {
    DoLDivEq(m,v);
    return v; 
  }

  template <class T, int N, StorageType S, IndexStyle I1, IndexStyle I2> 
  inline SmallVector<T,N,I1>& operator%=(
      SmallVector<T,N,I1>& v, const SmallMatrix<T,N,N,S,I2>& m)
  {
    DoRDivEq(m,v);
    return v; 
  }

  template <class T, int N, StorageType S, IndexStyle I1, IndexStyle I2> 
  inline SmallVector<CT,N,I1>& operator%=(
      SmallVector<CT,N,I1>& v, const SmallMatrix<T,N,N,S,I2>& m)
  { 
    DoRDivEq(m,v);
    return v; 
  }

  template <class T, int N, IndexStyle I1> 
  inline SmallVector<T,N,I1>& operator/=(
      SmallVector<T,N,I1>& v, const GenMatrix<T>& m)
  { 
    TMVAssert(m.rowsize() == N);
    TMVAssert(m.colsize() == N);
    m.LDivEq(v.View()); 
    return v; 
  }

  template <class T, int N, IndexStyle I1> 
  inline SmallVector<CT,N,I1>& operator/=(
      SmallVector<CT,N,I1>& v, const GenMatrix<T>& m)
  {
    TMVAssert(m.rowsize() == N);
    TMVAssert(m.colsize() == N);
    m.LDivEq(v.View()); 
    return v; 
  }

  template <class T, int N, IndexStyle I1> 
  inline SmallVector<T,N,I1>& operator%=(
      SmallVector<T,N,I1>& v, const GenMatrix<T>& m)
  {
    TMVAssert(m.rowsize() == N);
    TMVAssert(m.colsize() == N);
    m.RDivEq(v.View()); 
    return v; 
  }

  template <class T, int N, IndexStyle I1> 
  inline SmallVector<CT,N,I1>& operator%=(
      SmallVector<CT,N,I1>& v, const GenMatrix<T>& m)
  { 
    TMVAssert(m.rowsize() == N);
    TMVAssert(m.colsize() == N);
    m.RDivEq(v.View()); 
    return v; 
  }

  template <class T, int N, StorageType S, IndexStyle I2> 
  inline const VectorView<T>& operator/=(
      const VectorView<T>& v, const SmallMatrix<T,N,N,S,I2>& m)
  { 
    TMVAssert(v.size() == N);
    m.View().LDivEq(v); 
    return v; 
  }

  template <class T, int N, StorageType S, IndexStyle I2> 
  inline const VectorView<CT>& operator/=(
      const VectorView<CT>& v, const SmallMatrix<T,N,N,S,I2>& m)
  {
    TMVAssert(v.size() == N);
    m.View().LDivEq(v); 
    return v; 
  }

  template <class T, int N, StorageType S, IndexStyle I2> 
  inline const VectorView<T>& operator%=(
      const VectorView<T>& v, const SmallMatrix<T,N,N,S,I2>& m)
  {
    TMVAssert(v.size() == N);
    m.View().RDivEq(v); 
    return v; 
  }

  template <class T, int N, StorageType S, IndexStyle I2> 
  inline const VectorView<CT>& operator%=(
      const VectorView<CT>& v, const SmallMatrix<T,N,N,S,I2>& m)
  { 
    TMVAssert(v.size() == N);
    m.View().RDivEq(v); 
    return v; 
  }

  template <class T, class Tm, int N, StorageType S, IndexStyle I1, IndexStyle I2> 
  inline SmallVector<T,N,I1>& operator*=(
      SmallVector<T,N,I1>& v, const QuotXm_1<T,Tm,N,N,S,I2>& qxm)
  {
    DoRDivEq(qxm.GetM(),v);
    return v; 
  }

  template <class T, int N, StorageType S, IndexStyle I1, IndexStyle I2> 
  inline SmallVector<CT,N,I1>& operator*=(
      SmallVector<CT,N,I1>& v, const QuotXm_1<T,T,N,N,S,I2>& qxm)
  {
    DoRDivEq(qxm.GetM(),v);
    return v; 
  }

  template <class T, class Tm, int N, StorageType S, IndexStyle I1, IndexStyle I2> 
  inline SmallVector<T,N,I1>& operator*=(
      SmallVector<T,N,I1>& v, const QuotXm<T,Tm,N,N,S,I2>& qxm)
  {
    DoRDivEq(qxm.GetM(),v);
    MultXV<N>(qxm.GetX(),v.ptr());
    return v; 
  }

  template <class T, int N, StorageType S, IndexStyle I1, IndexStyle I2> 
  inline SmallVector<CT,N,I1>& operator*=(
      SmallVector<CT,N,I1>& v, const QuotXm<T,T,N,N,S,I2>& qxm)
  {
    DoRDivEq(qxm.GetM(),v);
    MultXV<N>(qxm.GetX(),v.ptr());
    return v; 
  }

  template <class T, class Tm, int N, StorageType S, IndexStyle I2> 
  inline const VectorView<T>& operator*=(
      const VectorView<T>& v, const QuotXm_1<T,Tm,N,N,S,I2>& qxm)
  {
    TMVAssert(v.size() == N);
    qxm.GetM().View().RDivEq(v); 
    return v; 
  }

  template <class T, int N, StorageType S, IndexStyle I2> 
  inline const VectorView<CT>& operator*=(
      const VectorView<CT>& v, const QuotXm_1<T,T,N,N,S,I2>& qxm)
  {
    TMVAssert(v.size() == N);
    qxm.GetM().View().RDivEq(v); 
    return v; 
  }

  template <class T, class Tm, int N, StorageType S, IndexStyle I2> 
  inline const VectorView<T>& operator*=(
      const VectorView<T>& v, const QuotXm<T,Tm,N,N,S,I2>& qxm)
  {
    TMVAssert(v.size() == N);
    qxm.GetM().View().RDivEq(v); 
    MultXV(qxm.GetX(),v.ptr());
    return v; 
  }

  template <class T, int N, StorageType S, IndexStyle I2> 
  inline const VectorView<CT>& operator*=(
      const VectorView<CT>& v, const QuotXm<T,T,N,N,S,I2>& qxm)
  {
    TMVAssert(v.size() == N);
    qxm.GetM().View().RDivEq(v); 
    MultXV(qxm.GetX(),v.ptr());
    return v; 
  }

  template <class T, class Tm, int N, IndexStyle I1> 
  inline SmallVector<T,N,I1>& operator*=(
      SmallVector<T,N,I1>& v, const QuotXM<T,Tm>& qxm)
  {
    TMVAssert(qxm.GetM().rowsize() == N);
    TMVAssert(qxm.GetM().colsize() == N);
    qxm.GetM().RDivEq(v.View()); 
    MultXV(qxm.GetX(),v.ptr());
    return v; 
  }

  template <class T, int N, IndexStyle I1> 
  inline SmallVector<CT,N,I1>& operator*=(
      SmallVector<CT,N,I1>& v, const QuotXM<T,T>& qxm)
  {
    TMVAssert(qxm.GetM().rowsize() == N);
    TMVAssert(qxm.GetM().colsize() == N);
    qxm.GetM().RDivEq(v.View()); 
    MultXV(qxm.GetX(),v.ptr());
    return v; 
  }

#define GENMATRIX1 SmallVector
#define GENMATRIX2 SmallMatrix
#define PRODXM1 ProdXv
#define PRODXM2 ProdXm
#define QUOTXM_1 QuotXm_1
#define QUOTXM QuotXm
#define QUOTMM_1 Quotvm_1
#define QUOTMM Quotvm
#define RQUOTMM_1 RQuotvm_1
#define RQUOTMM RQuotvm
#define X1 ,M,I1
#define X1b ,N,I1
#define X2 ,M,N,S,I2
#define X3 ,M,N,S,I1,I2
#define Y ,int M,int N,StorageType S, IndexStyle I1, IndexStyle I2
#define GETM1 .GetV()
#define GETM2 .GetM()
#include "tmv/TMV_AuxQuotMM.h"
#define PRODMM_1 Quotvm_1
#define PRODMM Quotvm
#define X3 ,M,N,S,I1,I2
#define Y ,int M,int N,StorageType S, IndexStyle I1, IndexStyle I2
#define GETM1 .GetV()
#define GETM2 .GetM()
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#define PRODMM_1 RQuotvm_1
#define PRODMM RQuotvm
#define X3 ,M,N,S,I1,I2
#define Y ,int M,int N,StorageType S, IndexStyle I1, IndexStyle I2
#define GETM1 .GetV()
#define GETM2 .GetM()
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTXM
#undef QUOTMM
#undef RQUOTMM

#define GENMATRIX1 SmallVector
#define GENMATRIX2 GenMatrix
#define PRODXM1 ProdXv
#define PRODXM2 ProdXM
#define QUOTXM QuotXM
#define QUOTMM QuotvM
#define RQUOTMM RQuotvM
#define X1 ,M,I1
#define X3 ,M,I1
#define Y ,int M, IndexStyle I1
#define GETM1 .GetV()
#define GETM2 .GetM()
#include "tmv/TMV_AuxQuotMM.h"
#define PRODMM QuotvM
#define X3 ,M,I1
#define Y ,int M, IndexStyle I1
#define GETM1 .GetV()
#define GETM2 .GetM()
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#define PRODMM RQuotvM
#define X3 ,M,I1
#define Y ,int M, IndexStyle I1
#define GETM1 .GetV()
#define GETM2 .GetM()
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTXM
#undef QUOTMM
#undef RQUOTMM

#define GENMATRIX1 GenVector
#define GENMATRIX2 SmallMatrix
#define PRODXM1 ProdXV
#define PRODXM2 ProdXm
#define QUOTXM_1 QuotXm_1
#define QUOTXM QuotXm
#define QUOTMM_1 QuotVm_1
#define QUOTMM QuotVm
#define RQUOTMM_1 RQuotVm_1
#define RQUOTMM RQuotVm
#define X2 ,M,N,S,I
#define X3 ,M,N,S,I
#define Y ,int M,int N,StorageType S,IndexStyle I
#define GETM1 .GetV()
#define GETM2 .GetM()
#include "tmv/TMV_AuxQuotMM.h"
#define PRODMM QuotVm
#define X3 ,M,N,S,I
#define Y ,int M,int N,StorageType S,IndexStyle I
#define GETM1 .GetV()
#define GETM2 .GetM()
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#define PRODMM RQuotVm
#define X3 ,M,N,S,I
#define Y ,int M,int N,StorageType S,IndexStyle I
#define GETM1 .GetV()
#define GETM2 .GetM()
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTXM
#undef QUOTMM
#undef RQUOTMM

  //
  // Matrix / % Matrix
  //

  template <class T, class T1, class T2, int M, int N, int K, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
  class Quotmm_1 :
    public SmallMatrixComposite<T,M,N>
  {
  public:
    inline Quotmm_1(
        const T DEBUGPARAM(_x), const SmallMatrix<T1,K,N,S1,I1>& _m1,
        const SmallMatrix<T2,K,M,S2,I2>& _m2) :
      m1(_m1), m2(_m2) { TMVAssert(_x == T(1)); }
    inline StorageType stor() const { return S1; }
    inline const SmallMatrix<T1,K,N,S1,I1>& GetM1() const { return m1; }
    inline const SmallMatrix<T2,K,M,S2,I2>& GetM2() const { return m2; }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,CStyle>& m0) const
    { TMVAssert(IsReal(T())); DoLDiv(m2,m1,m0); }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,CStyle>& m0) const
    { TMVAssert(IsReal(T())); DoLDiv(m2,m1,m0); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,CStyle>& m0) const
    { DoLDiv(m2,m1,m0); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,CStyle>& m0) const
    { DoLDiv(m2,m1,m0); }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,FortranStyle>& m0) const
    { TMVAssert(IsReal(T())); DoLDiv(m2,m1,m0); }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,FortranStyle>& m0) const
    { TMVAssert(IsReal(T())); DoLDiv(m2,m1,m0); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,FortranStyle>& m0) const
    { DoLDiv(m2,m1,m0); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,FortranStyle>& m0) const
    { DoLDiv(m2,m1,m0); }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == M && m0.rowsize() == N);
      TMVAssert(IsReal(T())); 
      m2.View().LDiv(m1.View(),m0);
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == M && m0.rowsize() == N);
      m2.View().LDiv(m1.View(),m0);
    }
  private:
    const SmallMatrix<T1,K,N,S1,I1>& m1;
    const SmallMatrix<T2,K,M,S2,I2>& m2;
  };

  template <class T, class T1, class T2, int M, int N, int K,StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
  class Quotmm :
    public SmallMatrixComposite<T,M,N>
  {
  public:
    inline Quotmm(const T _x, const SmallMatrix<T1,K,N,S1,I1>& _m1,
        const SmallMatrix<T2,K,M,S2,I2>& _m2) :
      x(_x), m1(_m1), m2(_m2) {}
    inline StorageType stor() const { return S1; }
    inline T GetX() const { return x; }
    inline const SmallMatrix<T1,K,N,S1,I1>& GetM1() const { return m1; }
    inline const SmallMatrix<T2,K,M,S2,I2>& GetM2() const { return m2; }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,CStyle>& m0) const
    { TMVAssert(IsReal(T())); DoLDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,CStyle>& m0) const
    { TMVAssert(IsReal(T())); DoLDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,CStyle>& m0) const
    { DoLDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,CStyle>& m0) const
    { DoLDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,FortranStyle>& m0) const
    { TMVAssert(IsReal(T())); DoLDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,FortranStyle>& m0) const
    { TMVAssert(IsReal(T())); DoLDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,FortranStyle>& m0) const
    { DoLDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,FortranStyle>& m0) const
    { DoLDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == M && m0.rowsize() == N);
      TMVAssert(IsReal(T())); 
      m2.View().LDiv(m1.View(),m0);
      MultXM(x,m0);
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == M && m0.rowsize() == N);
      m2.View().LDiv(m1.View(),m0);
      MultXM(x,m0);
    }
  private:
    const T x;
    const SmallMatrix<T1,K,N,S1,I1>& m1;
    const SmallMatrix<T2,K,M,S2,I2>& m2;
  };

  template <class T, class T1, class T2, int K, int M, StorageType S2, IndexStyle I2> 
  class QuotMm_1 :
    public MatrixComposite<T>
  {
  public:
    inline QuotMm_1(const T DEBUGPARAM(_x), const GenMatrix<T1>& _m1,
        const SmallMatrix<T2,K,M,S2,I2>& _m2) :
      m1(_m1), m2(_m2) { TMVAssert(_x==T(1)); TMVAssert(m1.colsize() == K); }
    inline size_t colsize() const { return M; }
    inline size_t rowsize() const { return m1.rowsize(); }
    inline StorageType stor() const { return S2; }
    inline const GenMatrix<T1>& GetM1() const { return m1; }
    inline const SmallMatrix<T2,K,M,S2,I2>& GetM2() const { return m2; }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == M && m0.rowsize() == rowsize());
      TMVAssert(IsReal(T())); 
      m2.View().LDiv(m1,m0);
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == M && m0.rowsize() == rowsize());
      m2.View().LDiv(m1,m0);
    }
  private:
    const GenMatrix<T1>& m1;
    const SmallMatrix<T2,K,M,S2,I2>& m2;
  };

  template <class T, class T1, class T2, int K, int M, StorageType S2, IndexStyle I2> 
  class QuotMm :
    public MatrixComposite<T>
  {
  public:
    inline QuotMm(const T _x, const GenMatrix<T1>& _m1,
        const SmallMatrix<T2,K,M,S2,I2>& _m2) :
      x(_x), m1(_m1), m2(_m2) { TMVAssert(m1.colsize() == K); }
    inline size_t colsize() const { return M; }
    inline size_t rowsize() const { return m1.rowsize(); }
    inline StorageType stor() const { return S2; }
    inline T GetX() const { return x; }
    inline const GenMatrix<T1>& GetM1() const { return m1; }
    inline const SmallMatrix<T2,K,M,S2,I2>& GetM2() const { return m2; }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == M && m0.rowsize() == rowsize());
      TMVAssert(IsReal(T())); 
      m2.View().LDiv(m1,m0);
      MultXM(x,m0);
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == M && m0.rowsize() == rowsize());
      m2.View().LDiv(m1,m0);
      MultXM(x,m0);
    }
  private:
    const T x;
    const GenMatrix<T1>& m1;
    const SmallMatrix<T2,K,M,S2,I2>& m2;
  };

  template <class T, class T1, class T2, int K, int N, StorageType S1, IndexStyle I1> 
  class QuotmM :
    public MatrixComposite<T>
  {
  public:
    inline QuotmM(const T _x, const SmallMatrix<T1,K,N,S1,I1>& _m1,
        const GenMatrix<T2>& _m2) :
      x(_x), m1(_m1), m2(_m2) {}
    inline size_t colsize() const { return m2.rowsize(); }
    inline size_t rowsize() const { return N; }
    inline StorageType stor() const { return S1; }
    inline T GetX() const { return x; }
    inline const SmallMatrix<T1,K,N,S1,I1>& GetM1() const { return m1; }
    inline const GenMatrix<T2>& GetM2() const { return m2; }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == colsize() && m0.rowsize() == N);
      TMVAssert(IsReal(T())); 
      m2.LDiv(m1.View(),m0);
      MultXM(x,m0);
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == colsize() && m0.rowsize() == N);
      m2.LDiv(m1.View(),m0);
      MultXM(x,m0);
    }
  private:
    const T x;
    const SmallMatrix<T1,K,N,S1,I1>& m1;
    const GenMatrix<T2>& m2;
  };

  template <class T, class T1, class T2, int M, int N, int K, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
  class RQuotmm_1 :
    public SmallMatrixComposite<T,M,N>
  {
  public:
    inline RQuotmm_1(
        const T DEBUGPARAM(_x), const SmallMatrix<T1,M,K,S1,I1>& _m1,
        const SmallMatrix<T2,N,K,S2,I2>& _m2) :
      m1(_m1), m2(_m2) { TMVAssert(_x == T(1)); }
    inline StorageType stor() const { return S1; }
    inline const SmallMatrix<T1,M,K,S1,I1>& GetM1() const { return m1; }
    inline const SmallMatrix<T2,N,K,S2,I2>& GetM2() const { return m2; }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,CStyle>& m0) const
    { TMVAssert(IsReal(T())); DoRDiv(m2,m1,m0); }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,CStyle>& m0) const
    { TMVAssert(IsReal(T())); DoRDiv(m2,m1,m0); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,CStyle>& m0) const
    { DoRDiv(m2,m1,m0); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,CStyle>& m0) const
    { DoRDiv(m2,m1,m0); }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,FortranStyle>& m0) const
    { TMVAssert(IsReal(T())); DoRDiv(m2,m1,m0); }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,FortranStyle>& m0) const
    { TMVAssert(IsReal(T())); DoRDiv(m2,m1,m0); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,FortranStyle>& m0) const
    { DoRDiv(m2,m1,m0); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,FortranStyle>& m0) const
    { DoRDiv(m2,m1,m0); }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == M && m0.rowsize() == N);
      TMVAssert(IsReal(T())); 
      m2.View().RDiv(m1.View(),m0);
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == M && m0.rowsize() == N);
      m2.View().RDiv(m1.View(),m0);
    }
  private:
    const SmallMatrix<T1,M,K,S1,I1>& m1;
    const SmallMatrix<T2,N,K,S2,I2>& m2;
  };

  template <class T, class T1, class T2, int M, int N, int K, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
  class RQuotmm :
    public SmallMatrixComposite<T,M,N>
  {
  public:
    inline RQuotmm(
        const T _x, const SmallMatrix<T1,M,K,S1,I1>& _m1,
        const SmallMatrix<T2,N,K,S2,I2>& _m2) :
      x(_x), m1(_m1), m2(_m2) {}
    inline StorageType stor() const { return S1; }
    inline const SmallMatrix<T1,M,K,S1,I1>& GetM1() const { return m1; }
    inline const SmallMatrix<T2,N,K,S2,I2>& GetM2() const { return m2; }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,CStyle>& m0) const
    { TMVAssert(IsReal(T())); DoRDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,CStyle>& m0) const
    { TMVAssert(IsReal(T())); DoRDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,CStyle>& m0) const
    { DoRDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,CStyle>& m0) const
    { DoRDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,ColMajor,FortranStyle>& m0) const
    { TMVAssert(IsReal(T())); DoRDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<RealType(T),M,N,RowMajor,FortranStyle>& m0) const
    { TMVAssert(IsReal(T())); DoRDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,ColMajor,FortranStyle>& m0) const
    { DoRDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
    inline void AssignTom(
        SmallMatrix<ComplexType(T),M,N,RowMajor,FortranStyle>& m0) const
    { DoRDiv(m2,m1,m0); MultXV<M*N>(x,m0.ptr()); }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == M && m0.rowsize() == N);
      TMVAssert(IsReal(T())); 
      m2.View().RDiv(m1.View(),m0);
      MultXM(x,m0);
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == M && m0.rowsize() == N);
      m2.View().RDiv(m1.View(),m0);
      MultXM(x,m0);
    }
  private:
    const T x;
    const SmallMatrix<T1,M,K,S1,I1>& m1;
    const SmallMatrix<T2,N,K,S2,I2>& m2;
  };

  template <class T, class T1, class T2, int N, int K, StorageType S2, IndexStyle I2> 
  class RQuotMm_1 :
    public MatrixComposite<T>
  {
  public:
    inline RQuotMm_1(const T DEBUGPARAM(_x), const GenMatrix<T1>& _m1,
        const SmallMatrix<T2,N,K,S2,I2>& _m2) :
      m1(_m1), m2(_m2) { TMVAssert(_x==T(1)); TMVAssert(m1.rowsize() == K); }
    inline size_t colsize() const { return m1.colsize(); }
    inline size_t rowsize() const { return N; }
    inline StorageType stor() const { return S2; }
    inline const GenMatrix<T1>& GetM1() const { return m1; }
    inline const SmallMatrix<T2,N,K,S2,I2>& GetM2() const { return m2; }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == colsize() && m0.rowsize() == N);
      TMVAssert(IsReal(T())); 
      m2.View().RDiv(m1,m0);
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == colsize() && m0.rowsize() == N);
      m2.View().RDiv(m1,m0);
    }
  private:
    const GenMatrix<T1>& m1;
    const SmallMatrix<T2,N,K,S2,I2>& m2;
  };

  template <class T, class T1, class T2, int N, int K, StorageType S2, IndexStyle I2> 
  class RQuotMm :
    public MatrixComposite<T>
  {
  public:
    inline RQuotMm(const T _x, const GenMatrix<T1>& _m1,
        const SmallMatrix<T2,N,K,S2,I2>& _m2) :
      x(_x), m1(_m1), m2(_m2) { TMVAssert(m1.rowsize() == K); }
    inline size_t colsize() const { return m1.colsize(); }
    inline size_t rowsize() const { return N; }
    inline StorageType stor() const { return S2; }
    inline T GetX() const { return x; }
    inline const GenMatrix<T1>& GetM1() const { return m1; }
    inline const SmallMatrix<T2,N,K,S2,I2>& GetM2() const { return m2; }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == colsize() && m0.rowsize() == N);
      TMVAssert(IsReal(T())); 
      m2.View().RDiv(m1,m0);
      MultXM(x,m0);
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == colsize() && m0.rowsize() == N);
      m2.View().RDiv(m1,m0);
      MultXM(x,m0);
    }
  private:
    const T x;
    const GenMatrix<T1>& m1;
    const SmallMatrix<T2,N,K,S2,I2>& m2;
  };

  template <class T, class T1, class T2, int M, int K, StorageType S1, IndexStyle I1> 
  class RQuotmM :
    public MatrixComposite<T>
  {
  public:
    inline RQuotmM(const T _x, const SmallMatrix<T1,M,K,S1,I1>& _m1,
        const GenMatrix<T2>& _m2) :
      x(_x), m1(_m1), m2(_m2) { TMVAssert(m2.rowsize() == K); }
    inline size_t colsize() const { return M; }
    inline size_t rowsize() const { return m2.colsize(); }
    inline StorageType stor() const { return S1; }
    inline T GetX() const { return x; }
    inline const SmallMatrix<T1,M,K,S1,I1>& GetM1() const { return m1; }
    inline const GenMatrix<T2>& GetM2() const { return m2; }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == M && m0.rowsize() == rowsize());
      TMVAssert(IsReal(T())); 
      m2.RDiv(m1.View(),m0);
      MultXM(x,m0);
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == M && m0.rowsize() == rowsize());
      m2.RDiv(m1.View(),m0);
      MultXM(x,m0);
    }
  private:
    const T x;
    const SmallMatrix<T1,M,K,S1,I1>& m1;
    const GenMatrix<T2>& m2;
  };

  template <class T, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
  inline SmallMatrix<T,M,N,S1,I1>& operator/=(
      SmallMatrix<T,M,N,S1,I1>& m1, const SmallMatrix<T,M,M,S2,I2>& m2)
  { 
    DoLDivEq(m2,m1); 
    return m1; 
  }

  template <class T, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
  inline SmallMatrix<CT,M,N,S1,I1>& operator/=(
      SmallMatrix<CT,M,N,S1,I1>& m1, const SmallMatrix<T,M,M,S2,I2>& m2)
  {
    DoLDivEq(m2,m1);
    return m1; 
  }

  template <class T, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
  inline SmallMatrix<T,M,N,S1,I1>& operator%=(
      SmallMatrix<T,M,N,S1,I1>& m1, const SmallMatrix<T,N,N,S2,I2>& m2)
  {
    DoRDivEq(m2,m1);
    return m1; 
  }

  template <class T, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
  inline SmallMatrix<CT,M,N,S1,I1>& operator%=(
      SmallMatrix<CT,M,N,S1,I1>& m1, const SmallMatrix<T,N,N,S2,I2>& m2)
  { 
    DoRDivEq(m2,m1);
    return m1; 
  }

  template <class T, int M, int N, StorageType S1, IndexStyle I1> 
  inline SmallMatrix<T,M,N,S1,I1>& operator/=(
      SmallMatrix<T,M,N,S1,I1>& m1, const GenMatrix<T>& m2)
  { 
    TMVAssert(m2.rowsize() == M && m2.colsize() == M);
    m2.LDivEq(m1.View()); 
    return m1; 
  }

  template <class T, int M, int N, StorageType S1, IndexStyle I1> 
  inline SmallMatrix<CT,M,N,S1,I1>& operator/=(
      SmallMatrix<CT,M,N,S1,I1>& m1, const GenMatrix<T>& m2)
  {
    TMVAssert(m2.rowsize() == M && m2.colsize() == M);
    m2.LDivEq(m1.View()); 
    return m1; 
  }

  template <class T, int M, int N, StorageType S1, IndexStyle I1> 
  inline SmallMatrix<T,M,N,S1,I1>& operator%=(
      SmallMatrix<T,M,N,S1,I1>& m1, const GenMatrix<T>& m2)
  {
    TMVAssert(m2.rowsize() == N && m2.colsize() == N);
    m2.RDivEq(m1.View()); 
    return m1; 
  }

  template <class T, int M, int N, StorageType S1, IndexStyle I1> 
  inline SmallMatrix<CT,M,N,S1,I1>& operator%=(
      SmallMatrix<CT,M,N,S1,I1>& m1, const GenMatrix<T>& m2)
  { 
    TMVAssert(m2.rowsize() == N && m2.colsize() == N);
    m2.RDivEq(m1.View()); 
    return m1; 
  }

  template <class T, int M, StorageType S2, IndexStyle I2> 
  inline const MatrixView<T>& operator/=(
      const MatrixView<T>& m1, const SmallMatrix<T,M,M,S2,I2>& m2)
  { 
    TMVAssert(m1.colsize() == M);
    m2.View().LDivEq(m1); 
    return m1; 
  }

  template <class T, int M, StorageType S2, IndexStyle I2> 
  inline const MatrixView<CT>& operator/=(
      const MatrixView<CT>& m1, const SmallMatrix<T,M,M,S2,I2>& m2)
  {
    TMVAssert(m1.colsize() == M);
    m2.View().LDivEq(m1); 
    return m1; 
  }

  template <class T, int N, StorageType S2, IndexStyle I2> 
  inline const MatrixView<T>& operator%=(
      const MatrixView<T>& m1, const SmallMatrix<T,N,N,S2,I2>& m2)
  {
    TMVAssert(m1.rowsize() == N);
    m2.View().RDivEq(m1); 
    return m1; 
  }

  template <class T, int N, StorageType S2, IndexStyle I2> 
  inline const MatrixView<CT>& operator%=(
      const MatrixView<CT>& m1, const SmallMatrix<T,N,N,S2,I2>& m2)
  { 
    TMVAssert(m1.rowsize() == N);
    m2.View().RDivEq(m1); 
    return m1; 
  }

  template <class T, class Tm, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
  inline SmallMatrix<T,M,N,S1,I1>& operator*=(
      SmallMatrix<T,M,N,S1,I1>& m1, const QuotXm_1<T,Tm,N,N,S2,I2>& qxm)
  {
    DoRDivEq(qxm.GetM(),m1); 
    return m1; 
  }

  template <class T, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
  inline SmallMatrix<CT,M,N,S1,I1>& operator*=(
      SmallMatrix<CT,M,N,S1,I1>& m1, const QuotXm_1<T,T,N,N,S2,I2>& qxm)
  {
    DoRDivEq(qxm.GetM(),m1);
    return m1; 
  }

  template <class T, class Tm, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
  inline SmallMatrix<T,M,N,S1,I1>& operator*=(
      SmallMatrix<T,M,N,S1,I1>& m1, const QuotXm<T,Tm,N,N,S2,I2>& qxm)
  {
    DoRDivEq(qxm.GetM(),m1);
    MultXV<M*N>(qxm.GetX(),m1.ptr());
    return m1; 
  }

  template <class T, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
  inline SmallMatrix<CT,M,N,S1,I1>& operator*=(
      SmallMatrix<CT,M,N,S1,I1>& m1, const QuotXm<T,T,N,N,S2,I2>& qxm)
  {
    DoRDivEq(qxm.GetM(),m1);
    MultXV<M*N>(qxm.GetX(),m1.ptr());
    return m1; 
  }

  template <class T, class Tm, int N, StorageType S2, IndexStyle I2> 
  inline const MatrixView<T>& operator*=(
      const MatrixView<T>& m1, const QuotXm_1<T,Tm,N,N,S2,I2>& qxm)
  {
    TMVAssert(m1.rowsize() == N);
    qxm.GetM().View().RDivEq(m1); 
    return m1; 
  }

  template <class T, int N, StorageType S2, IndexStyle I2> 
  inline const MatrixView<CT>& operator*=(
      const MatrixView<CT>& m1, const QuotXm_1<T,T,N,N,S2,I2>& qxm)
  {
    TMVAssert(m1.rowsize() == N);
    qxm.GetM().View().RDivEq(m1); 
    return m1; 
  }

  template <class T, class Tm, int N, StorageType S2, IndexStyle I2> 
  inline const MatrixView<T>& operator*=(
      const MatrixView<T>& m1, const QuotXm<T,Tm,N,N,S2,I2>& qxm)
  {
    TMVAssert(m1.rowsize() == N);
    qxm.GetM().View().RDivEq(m1); 
    m1 *= qxm.GetX();
    return m1; 
  }

  template <class T, int N, StorageType S2, IndexStyle I2> 
  inline const MatrixView<CT>& operator*=(
      const MatrixView<CT>& m1, const QuotXm<T,T,N,N,S2,I2>& qxm)
  {
    TMVAssert(m1.rowsize() == N);
    qxm.GetM().View().RDivEq(m1); 
    m1 *= qxm.GetX();
    return m1; 
  }

  template <class T, class Tm, int M, int N, StorageType S1, IndexStyle I1> 
  inline SmallMatrix<T,M,N,S1,I1>& operator*=(
      SmallMatrix<T,M,N,S1,I1>& m1, const QuotXM<T,Tm>& qxm)
  {
    TMVAssert(qxm.GetM().rowsize() == N);
    TMVAssert(qxm.GetM().colsize() == N);
    qxm.GetM().RDivEq(m1.View()); 
    m1 *= qxm.GetX();
    return m1; 
  }

  template <class T, int M, int N, StorageType S1, IndexStyle I1> 
  inline SmallMatrix<CT,M,N,S1,I1>& operator*=(
      SmallMatrix<CT,M,N,S1,I1>& m1, const QuotXM<T,T>& qxm)
  {
    TMVAssert(qxm.GetM().rowsize() == N);
    TMVAssert(qxm.GetM().colsize() == N);
    qxm.GetM().RDivEq(m1.View()); 
    m1 *= qxm.GetX();
    return m1; 
  }

#define GENMATRIX1 SmallMatrix
#define GENMATRIX2 SmallMatrix
#define PRODXM1 ProdXm
#define PRODXM2 ProdXm
#define QUOTXM_1 QuotXm_1
#define QUOTXM QuotXm
#define QUOTMM_1 Quotmm_1
#define QUOTMM Quotmm
#define RQUOTMM_1 RQuotmm_1
#define RQUOTMM RQuotmm
#define X1 ,K,N,S1,I1
#define X1b ,M,K,S1,I1
#define X2 ,K,M,S2,I2
#define X2b ,N,K,S2,I2
#define X3 ,M,N,K,S1,S2,I1,I2
#define Y ,int M,int N,int K,StorageType S1,StorageType S2,IndexStyle I1,IndexStyle I2
#include "tmv/TMV_AuxQuotMM.h"
#define PRODMM_1 Quotmm_1
#define PRODMM Quotmm
#define X3 ,M,N,K,S1,S2,I1,I2
#define Y ,int M,int N,int K,StorageType S1,StorageType S2,IndexStyle I1,IndexStyle I2
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#define PRODMM_1 RQuotmm_1
#define PRODMM RQuotmm
#define X3 ,M,N,K,S1,S2,I1,I2
#define Y ,int M,int N,int K,StorageType S1,StorageType S2,IndexStyle I1,IndexStyle I2
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTXM
#undef QUOTMM
#undef RQUOTMM

#define GENMATRIX1 SmallMatrix
#define GENMATRIX2 GenMatrix
#define PRODXM1 ProdXm
#define PRODXM2 ProdXM
#define QUOTXM QuotXM
#define QUOTMM QuotmM
#define RQUOTMM RQuotmM
#define X1 ,K,N,S1,I1
#define X1b ,M,K,S1,I1
#define X3 ,K,N,S1,I1
#define X3b ,M,K,S1,I1
#define Y ,int K, int N, StorageType S1, IndexStyle I1
#define Yb ,int M, int K, StorageType S1, IndexStyle I1
#include "tmv/TMV_AuxQuotMM.h"
#define PRODMM QuotmM
#define X3 ,K,N,S1,I1
#define Y ,int K, int N, StorageType S1, IndexStyle I1
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#define PRODMM RQuotmM
#define X3 ,M,K,S1,I1
#define Y ,int M, int K, StorageType S1, IndexStyle I1
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTXM
#undef QUOTMM
#undef RQUOTMM

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 SmallMatrix
#define PRODXM1 ProdXM
#define PRODXM2 ProdXm
#define QUOTXM_1 QuotXm_1
#define QUOTXM QuotXm
#define QUOTMM_1 QuotMm_1
#define QUOTMM QuotMm
#define RQUOTMM_1 RQuotMm_1
#define RQUOTMM RQuotMm
#define X2 ,K,M,S2,I2
#define X2b ,N,K,S2,I2
#define X3 ,K,M,S2,I2
#define X3b ,N,K,S2,I2
#define Y ,int K,int M,StorageType S2, IndexStyle I2
#define Yb ,int N,int K,StorageType S2, IndexStyle I2
#include "tmv/TMV_AuxQuotMM.h"
#define PRODMM QuotMm
#define X3 ,K,M,S2,I2
#define Y ,int K,int M,StorageType S2, IndexStyle I2
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#define PRODMM RQuotMm
#define X3 ,N,K,S2,I2
#define Y ,int N,int K,StorageType S2, IndexStyle I2
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTXM
#undef QUOTMM
#undef RQUOTMM

} // namespace tmv

#undef CT
#undef CCT
#undef VCT

#endif
