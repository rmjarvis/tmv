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


#ifndef TMV_MatrixArith_H
#define TMV_MatrixArith_H

#include "TMV_VectorArithFunc.h"
#include "TMV_MatrixArithFunc.h"

#define CT std::complex<T>
#define CCT ConjRef<std::complex<T> >
#define VCT VarConjRef<std::complex<T> >

namespace tmv {

  template <class T, class Tv> class ProdXV;

  template <class T, StorageType S, IndexStyle I, class Tx> 
    inline Matrix<T,S>& operator+=(Matrix<T,S,I>& m, const Tx& x) 
    { m.View() += x; return m; }

  template <class T, StorageType S, IndexStyle I, class Tx> 
    inline Matrix<T,S>& operator-=(Matrix<T,S,I>& m, const Tx& x) 
    { m.View() -= x; return m; }

  template <class T, StorageType S, IndexStyle I, class Tx> 
    inline Matrix<T,S>& operator*=(Matrix<T,S,I>& m, const Tx& x) 
    { m.View() *= x; return m; }

  template <class T, StorageType S, IndexStyle I, class Tx> 
    inline Matrix<T,S>& operator/=(Matrix<T,S,I>& m, const Tx& x) 
    { m.View() /= x; return m; }

  template <class T, StorageType S, IndexStyle I, class Tx> 
    inline Matrix<T,S>& operator%=(Matrix<T,S,I>& m, const Tx& x) 
    { m.View() %= x; return m; }

  //
  // Scalar * Matrix
  //

  template <class T, class Tm> class ProdXM : 
    public MatrixComposite<T> 
  {
    // x*m
    public:
      inline ProdXM(const T _x, const GenMatrix<Tm>& _m) : x(_x), m(_m) {}
      inline size_t colsize() const { return m.colsize(); }
      inline size_t rowsize() const { return m.rowsize(); }
      inline StorageType stor() const { return BaseStorOf(m); }
      inline T GetX() const { return x; }
      inline const GenMatrix<Tm>& GetM() const { return m; }
      inline void AssignToM(const MatrixView<RealType(T)>& m0) const
      {
	TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
	TMVAssert(IsReal(T()));
	MultXM(x,m0=m);
      }
      inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
	MultXM(x,m0=m);
      }
    private:
      const T x;
      const GenMatrix<Tm>& m;
  };

  // m*=x
  template <class T> inline const MatrixView<T>& operator*=(
      const MatrixView<T>& m, T x) 
  { MultXM(x,m); return m; }

  template <class T> inline const MatrixView<CT>& operator*=(
      const MatrixView<CT>& m, T x) 
  { MultXM(T(x),m); return m; }

  template <class T> inline const MatrixView<CT>& operator*=(
      const MatrixView<CT>& m, CCT x) 
  { MultXM(CT(x),m); return m; }

  template <class T> inline const MatrixView<CT>& operator*=(
      const MatrixView<CT>& m, VCT x) 
  { MultXM(CT(x),m); return m; }

  // m/=x
  template <class T> inline const MatrixView<T>& operator/=(
      const MatrixView<T>& m, T x) 
  { MultXM(T(1)/x,m); return m; }

  template <class T> inline const MatrixView<CT>& operator/=(
      const MatrixView<CT>& m, T x) 
  { MultXM(T(1)/x,m); return m; }

  template <class T> inline const MatrixView<CT>& operator/=(
      const MatrixView<CT>& m, CCT x) 
  { MultXM(T(1)/CT(x),m); return m; }

  template <class T> inline const MatrixView<CT>& operator/=(
      const MatrixView<CT>& m, VCT x) 
  { MultXM(T(1)/CT(x),m); return m; }

#define GENMATRIX GenMatrix
#define PRODXM ProdXM
#include "TMV_AuxProdXM.h"
  // Defines things like -m, x*m, m*x, x*(x*m), etc.
#undef GENMATRIX
#undef PRODXM

  //
  // Matrix + Scalar
  //

  template <class T, class Tm> class SumMX : 
    public MatrixComposite<T>
  {
    // x1*m + x2
    public:
      inline SumMX(T _x1, const GenMatrix<Tm>& _m, T _x2) :
	x1(_x1), m(_m), x2(_x2)
      { TMVAssert(m.IsSquare()); }
      inline size_t colsize() const { return m.colsize(); }
      inline size_t rowsize() const { return m.rowsize(); }
      inline StorageType stor() const { return BaseStorOf(m); }
      inline T GetX1() const { return x1; }
      inline const GenMatrix<Tm>& GetM() const { return m; }
      inline T GetX2() const { return x2; }
      inline void AssignToM(const MatrixView<RealType(T)>& m0) const
      { 
	TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
	TMVAssert(IsReal(T()));
	m0 = x1*m;
	m0.diag().AddToAll(REAL(x2));
      }
      inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
      { 
	TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
	m0 = x1*m;
	m0.diag().AddToAll(ComplexType(T)(x2));
      }
    private:
      const T x1;
      const GenMatrix<Tm>& m;
      const T x2;
  };

  // m+=x
  template <class T> inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m, T x) 
  {
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(x); 
    return m; 
  }

  template <class T> inline const MatrixView<CT>& operator+=(
      const MatrixView<CT>& m, T x) 
  {
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(CT(x)); 
    return m; 
  }

  template <class T> inline const MatrixView<CT>& operator+=(
      const MatrixView<CT>& m, CCT x) 
  {
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(CT(x)); 
    return m; 
  }

  template <class T> inline const MatrixView<CT>& operator+=(
      const MatrixView<CT>& m, VCT x) 
  {
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(CT(x)); 
    return m;
  }

  // m-=x
  template <class T> inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m, T x) 
  {
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(-x); 
    return m; 
  }

  template <class T> inline const MatrixView<CT>& operator-=(
      const MatrixView<CT>& m, T x) 
  {
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(CT(-x)); 
    return m; 
  }

  template <class T> inline const MatrixView<CT>& operator-=(
      const MatrixView<CT>& m, CCT x) 
  {
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(-CT(x)); 
    return m; 
  }

  template <class T> inline const MatrixView<CT>& operator-=(
      const MatrixView<CT>& m, VCT x) 
  {
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(-CT(x)); 
    return m; 
  }

#define GENMATRIX GenMatrix
#define PRODXM ProdXM
#define SUMMX SumMX
#include "TMV_AuxSumMX.h"
  // Defines things like m+x, x+m, x-m, m-x, x+x*m, x*(x+m), etc.
#undef SUMMX
#undef GENMATRIX
#undef PRODXM

  //
  // Vector ^ Vector (OuterProduct)
  //

  template <class T, class T1, class T2> class OProdVV : 
    public MatrixComposite<T>
  {
    public:
      inline OProdVV(const T _x, const GenVector<T1>& _v1,
	  const GenVector<T2>& _v2) :
	x(_x), v1(_v1), v2(_v2) {}
      inline size_t colsize() const { return v1.size(); }
      inline size_t rowsize() const { return v2.size(); }
      inline StorageType stor() const { return RowMajor; } // arbitrary
      inline T GetX() const { return x; }
      inline const GenVector<T1>& GetV1() const { return v1; }
      inline const GenVector<T2>& GetV2() const { return v2; }
      inline void AssignToM(const MatrixView<RealType(T)>& m0) const
      {
	TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
	TMVAssert(IsReal(T()));
	Rank1Update<false>(x, v1, v2, m0);
      }
      inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
	Rank1Update<false>(x, v1, v2, m0);
      }

    private:
      T x;
      const GenVector<T1>& v1;
      const GenVector<T2>& v2;
  };

  // m+=(x*v^v)
  template <class T, class T1, class T2> inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m0, const OProdVV<T,T1,T2>& opvv)
  { 
    TMVAssert(m0.colsize() == opvv.colsize());
    TMVAssert(m0.rowsize() == opvv.rowsize());
    Rank1Update<true>(opvv.GetX(), opvv.GetV1(), opvv.GetV2(), m0);
    return m0;
  }

  template <class T> inline const MatrixView<CT>& operator+=(
      const MatrixView<CT>& m0, const OProdVV<T,T,T>& opvv)
  { 
    TMVAssert(m0.colsize() == opvv.colsize());
    TMVAssert(m0.rowsize() == opvv.rowsize());
    Rank1Update<true>(opvv.GetX(), opvv.GetV1(), opvv.GetV2(), m0);
    return m0;
  }

  // m-=(x*v^v)
  template <class T, class T1, class T2> inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m0, const OProdVV<T,T1,T2>& opvv)
  { 
    TMVAssert(m0.colsize() == opvv.colsize());
    TMVAssert(m0.rowsize() == opvv.rowsize());
    Rank1Update<true>(-opvv.GetX(), opvv.GetV1(), opvv.GetV2(), m0);
    return m0;
  }

  template <class T> inline const MatrixView<CT>& operator-=(
      const MatrixView<CT>& m0, const OProdVV<T,T,T>& opvv)
  { 
    TMVAssert(m0.colsize() == opvv.colsize());
    TMVAssert(m0.rowsize() == opvv.rowsize());
    Rank1Update<true>(-opvv.GetX(), opvv.GetV1(), opvv.GetV2(), m0);
    return m0;
  }

#define OPRODVV OProdVV
#define GENVECTOR1 GenVector
#define GENVECTOR2 GenVector
#define PRODXV1 ProdXV
#define PRODXV2 ProdXV
#include "TMV_AuxOProdVV.h"
#undef OPRODVV
#undef GENVECTOR1
#undef GENVECTOR2
#undef PRODXV1
#undef PRODXV2


  //
  // Matrix + Matrix
  //

  template <class T, class T1, class T2> class SumMM : 
    public MatrixComposite<T> 
  {
    public:
      inline SumMM(const T _x1, const GenMatrix<T1>& _m1, 
	  const T _x2, const GenMatrix<T2>& _m2) :
	x1(_x1), m1(_m1), x2(_x2), m2(_m2)
      { 
	TMVAssert(m1.colsize() == m2.colsize());
	TMVAssert(m1.rowsize() == m2.rowsize()); 
      }
      inline size_t colsize() const { return m1.colsize(); }
      inline size_t rowsize() const { return m1.rowsize(); }
      inline StorageType stor() const { return BaseStorOf(m1); }
      inline T GetX1() const { return x1; }
      inline const GenMatrix<T1>& GetM1() const { return m1; }
      inline T GetX2() const { return x2; }
      inline const GenMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToM(const MatrixView<RealType(T)>& m0) const
      { 
	TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
	TMVAssert(IsReal(T()));
	AddMM(x1,m1,x2,m2,m0);
      }
      inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
      { 
	TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
	AddMM(x1,m1,x2,m2,m0);
      }
    private:
      const T x1;
      const GenMatrix<T1>& m1;
      const T x2;
      const GenMatrix<T2>& m2;
  };

  template <class T> inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m1, const GenMatrix<T>& m2) 
  { 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    AddMM(T(1),m2,m1); 
    return m1; 
  }

  template <class T> inline const MatrixView<CT>& operator+=(
      const MatrixView<CT>& m1, const GenMatrix<T>& m2) 
  {
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    AddMM(T(1),m2,m1); 
    return m1; 
  }

  template <class T> inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m1, const GenMatrix<T>& m2) 
  { 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    AddMM(T(-1),m2,m1); 
    return m1; 
  }

  template <class T> inline const MatrixView<CT>& operator-=(
      const MatrixView<CT>& m1, const GenMatrix<T>& m2) 
  { 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    AddMM(T(-1),m2,m1); 
    return m1; 
  }

  template <class T, class T2> inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m, const ProdXM<T,T2>& pxm)
  {
    TMVAssert(m.colsize() == pxm.colsize());
    TMVAssert(m.rowsize() == pxm.rowsize());
    AddMM(pxm.GetX(),pxm.GetM(),m);
    return m;
  }

  template <class T> inline const MatrixView<CT>& operator+=(
      const MatrixView<CT>& m, const ProdXM<T,T>& pxm)
  {
    TMVAssert(m.colsize() == pxm.colsize());
    TMVAssert(m.rowsize() == pxm.rowsize());
    AddMM(pxm.GetX(),pxm.GetM(),m);
    return m;
  }

  template <class T, class T2> inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m, const ProdXM<T,T2>& pxm)
  {
    TMVAssert(m.colsize() == pxm.colsize());
    TMVAssert(m.rowsize() == pxm.rowsize());
    AddMM(-pxm.GetX(),pxm.GetM(),m);
    return m;
  }

  template <class T> inline const MatrixView<CT>& operator-=(
      const MatrixView<CT>& m, const ProdXM<T,T>& pxm)
  {
    TMVAssert(m.colsize() == pxm.colsize());
    TMVAssert(m.rowsize() == pxm.rowsize());
    AddMM(-pxm.GetX(),pxm.GetM(),m);
    return m;
  }

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 GenMatrix
#define PRODXM1 ProdXM
#define PRODXM2 ProdXM
#define SUMMM SumMM
#include "TMV_AuxSumMM.h"
  // Defines things like m+m, m-m, (x*m)-(x*m), etc.
#include "TMV_AuxSumMMa.h"
  // Defines things like -(m+m), x*(m+m), (m+m)/x, etc.
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef SUMMM

  //
  // Matrix * Vector
  //

  template <class T, class T1, class T2> class ProdMV :
    public VectorComposite<T>
  {
    public:
      inline ProdMV(const T _x, const GenMatrix<T1>& _m,
	  const GenVector<T2>& _v) : x(_x), m(_m), v(_v)
      { TMVAssert(v.size()==m.rowsize()); }
      inline size_t size() const { return m.colsize(); }
      inline T GetX() const { return x; }
      inline const GenMatrix<T1>& GetM() const { return m; }
      inline const GenVector<T2>& GetV() const { return v; }
      inline void AssignToV(const VectorView<RealType(T)>& v0) const
      {
	TMVAssert(v0.size() == size());
	TMVAssert(IsReal(T()));
	MultMV<false>(x,m,v,v0);
      }
      inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
      {
	TMVAssert(v0.size() == size());
	MultMV<false>(x,m,v,v0);
      }
    private:
      const T x;
      const ConstMatrixView<T1> m;
      const GenVector<T2>& v;
  };

  template <class T> inline const VectorView<T>& operator*=(
      const VectorView<T>& v, const GenMatrix<T>& m)
  { 
    TMVAssert(v.size() == m.colsize());
    MultMV<false>(T(1),m.Transpose(),v,v);
    return v;
  }

  template <class T> inline const VectorView<CT>& operator*=(
      const VectorView<CT>& v, const GenMatrix<T>& m)
  {
    TMVAssert(v.size() == m.colsize());
    MultMV<false>(T(1),m.Transpose(),v,v);
    return v;
  }

  template <class T, class T2> inline const VectorView<T>& operator*=(
      const VectorView<T>& v, const ProdXM<T,T2>& pxm)
  { 
    TMVAssert(v.size() == pxm.colsize());
    MultMV<false>(pxm.GetX(),pxm.GetM().Transpose(),v,v);
    return v;
  }

  template <class T> inline const VectorView<CT>& operator*=(
      const VectorView<CT>& v, const ProdXM<T,T>& pxm)
  {
    TMVAssert(v.size() == pxm.colsize());
    MultMV<false>(pxm.GetX(),pxm.GetM().Transpose(),v,v);
    return v;
  }

  template <class T, class T1, class T2> inline const VectorView<T>& operator+=(
      const VectorView<T>& v, const ProdMV<T,T1,T2>& pmv)
  { 
    TMVAssert(v.size() == pmv.size());
    MultMV<true>(pmv.GetX(),pmv.GetM(),pmv.GetV(),v);
    return v;
  }

  template <class T> inline const VectorView<CT>& operator+=(
      const VectorView<CT>& v, const ProdMV<T,T,T>& pmv)
  { 
    TMVAssert(v.size() == pmv.size());
    MultMV<true>(pmv.GetX(),pmv.GetM(),pmv.GetV(),v);
    return v;
  }

  template <class T, class T1, class T2> inline const VectorView<T>& operator-=(
      const VectorView<T>& v, const ProdMV<T,T1,T2>& pmv)
  { 
    TMVAssert(v.size() == pmv.size());
    MultMV<true>(-pmv.GetX(), pmv.GetM(), pmv.GetV(), v);
    return v;
  }

  template <class T> inline const VectorView<CT>& operator-=(
      const VectorView<CT>& v, const ProdMV<T,T,T>& pmv)
  { 
    TMVAssert(v.size() == pmv.size());
    MultMV<true>(-pmv.GetX(),pmv.GetM(),pmv.GetV(),v);
    return v;
  }

#define GENMATRIX GenMatrix
#define GENVECTOR GenVector
#define PRODXV ProdXV
#define PRODXM ProdXM
#define PRODMV ProdMV
#include "TMV_AuxProdMV.h"
  // Defines things like v*m, m*v, (x*m)*(x*v), etc.
#include "TMV_AuxProdMVa.h"
  // Defines things like -(m*v), x*(m*v), (m*v)/x, etc.
#undef GENMATRIX
#undef GENVECTOR
#undef PRODXV
#undef PRODXM
#undef PRODMV

  //
  // Matrix * Matrix
  //

  template <class T, class T1, class T2> class ProdMM : 
    public MatrixComposite<T>
  {
    public:
      inline ProdMM(const T _x, const GenMatrix<T1>& _m1, 
	  const GenMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert(m1.rowsize() == m2.colsize()) ; }
      inline size_t colsize() const { return m1.colsize(); }
      inline size_t rowsize() const { return m2.rowsize(); }
      inline StorageType stor() const { return BaseStorOf(m1); }
      inline T GetX() const { return x; }
      inline const GenMatrix<T1>& GetM1() const { return m1; }
      inline const GenMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToM(const MatrixView<RealType(T)>& m0) const
      {
	TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
	TMVAssert(IsReal(T()));
	MultMM<false>(x, m1, m2, m0);
      }
      inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
	MultMM<false>(x, m1, m2, m0);
      }

    private:
      const T x;
      const GenMatrix<T1>& m1;
      const GenMatrix<T2>& m2;
  };

  template <class T> inline const MatrixView<T>& operator*=(
      const MatrixView<T>& m1, const GenMatrix<T>& m2)
  { 
    TMVAssert(m2.colsize()==m2.rowsize());
    TMVAssert(m1.rowsize()==m2.rowsize());
    MultMM<false>(T(1),m1,m2,m1); 
    return m1; 
  }

  template <class T> inline const MatrixView<CT>& operator*=(
      const MatrixView<CT>& m1, const GenMatrix<T>& m2)
  { 
    TMVAssert(m2.colsize()==m2.rowsize());
    TMVAssert(m1.rowsize()==m2.rowsize());
    MultMM<false>(T(1),m1,m2,m1); 
    return m1; 
  }

  template <class T, class T2> inline const MatrixView<T>& operator*=(
      const MatrixView<T>& m1, const ProdXM<T,T2>& pxm)
  { 
    TMVAssert(pxm.colsize()==pxm.rowsize());
    TMVAssert(m1.rowsize()==pxm.rowsize());
    MultMM<false>(pxm.GetX(),m1,pxm.GetM(),m1); 
    return m1; 
  }

  template <class T> inline const MatrixView<CT>& operator*=(
      const MatrixView<CT>& m1, const ProdXM<T,T>& pxm)
  { 
    TMVAssert(pxm.colsize()==pxm.rowsize());
    TMVAssert(m1.rowsize()==pxm.rowsize());
    MultMM<false>(pxm.GetX(),m1,pxm.GetM(),m1); 
    return m1; 
  }

  template <class T, class T2, class T3> inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m, const ProdMM<T,T2,T3>& pmm)
  { 
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
    return m; 
  }

  template <class T> inline const MatrixView<CT>& operator+=(
      const MatrixView<CT>& m, const ProdMM<T,T,T>& pmm)
  { 
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
    return m; 
  }

  template <class T, class T2, class T3> inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m, const ProdMM<T,T2,T3>& pmm)
  { 
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
    return m; 
  }

  template <class T> inline const MatrixView<CT>& operator-=(
      const MatrixView<CT>& m, const ProdMM<T,T,T>& pmm)
  {
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
    return m; 
  }

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 GenMatrix
#define PRODXM1 ProdXM
#define PRODXM2 ProdXM
#define PRODMM ProdMM
#include "TMV_AuxProdMM.h"
  // Defines things like m*m, m*(x*m), etc.
#include "TMV_AuxProdMMa.h"
  // Defines things like -(m*m), x*(m*m), etc.
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef PRODMM


  //
  // Scalar / Matrix
  //

  template <class T, class Tm> class QuotXM : 
    public MatrixComposite<T> 
  {
    public:
      inline QuotXM(const T _x, const GenMatrix<Tm>& _m) : x(_x), m(_m) {}
      inline size_t colsize() const { return m.rowsize(); }
      inline size_t rowsize() const { return m.colsize(); }
      inline StorageType stor() const { return ColMajor; }
      inline T GetX() const { return x; }
      inline const GenMatrix<Tm>& GetM() const { return m; }
      inline void AssignToM(const MatrixView<RealType(T)>& m0) const
      {
	TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
	TMVAssert(IsReal(T()));
	m.Inverse(m0);
	MultXM(x,m0);
      }
      inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
	m.Inverse(m0);
	MultXM(x,m0);
      }
    private:
      const T x;
      const GenMatrix<Tm>& m;
  };

#define GENMATRIX GenMatrix
#define PRODXM ProdXM
#define QUOTXM QuotXM
#include "TMV_AuxQuotXM.h"
#undef GENMATRIX
#undef PRODXM
#undef QUOTXM


  //
  // Vector / % Matrix 
  // v/m is the solution (x) of mx = v
  // ie. / is really division from the left: x = m^-1 v
  // Use % if you want division from the right (v m^-1)
  //

  template <class T, class T1, class T2> class QuotVM : 
    public VectorComposite<T>
  {
    public:
      inline QuotVM(const T _x, const GenVector<T1>& _v,
	  const GenMatrix<T2>& _m) :
	x(_x), v(_v), m(_m)
      { TMVAssert(v.size()==m.colsize()); }
      inline size_t size() const { return m.rowsize(); }
      inline T GetX() const { return x; }
      inline const GenVector<T1>& GetV() const { return v; }
      inline const GenMatrix<T2>& GetM() const { return m; }
      inline void AssignToV(const VectorView<RealType(T)>& v0) const
      {
	TMVAssert(v0.size() == size());
	TMVAssert(IsReal(T()));
	m.LDiv(v,v0);
	MultXV(x,v0);
      }
      inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
      {
	TMVAssert(v0.size() == size());
	m.LDiv(v,v0);
	MultXV(x,v0);
      }
    private:
      const T x;
      const GenVector<T1>& v;
      const GenMatrix<T2>& m;
  };

  template <class T, class T1, class T2> class RQuotVM : 
    public VectorComposite<T>
  {
    public:
      inline RQuotVM(const T _x, const GenVector<T1>& _v,
	  const GenMatrix<T2>& _m) :
	x(_x), v(_v), m(_m)
      { TMVAssert(v.size()==m.rowsize()); }
      inline size_t size() const { return m.colsize(); }
      inline T GetX() const { return x; }
      inline const GenVector<T1>& GetV() const { return v; }
      inline const GenMatrix<T2>& GetM() const { return m; }
      inline void AssignToV(const VectorView<RealType(T)>& v0) const
      {
	TMVAssert(v0.size() == size());
	TMVAssert(IsReal(T()));
	m.RDiv(v,v0);
	MultXV(x,v0);
      }
      inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
      {
	TMVAssert(v0.size() == size());
	m.RDiv(v,v0);
	MultXV(x,v0);
      }
    private:
      const T x;
      const GenVector<T1>& v;
      const GenMatrix<T2>& m;
  };

  template <class T> inline const VectorView<T>& operator/=(
      const VectorView<T>& v, const GenMatrix<T>& m)
  { 
    TMVAssert(m.IsSquare());
    TMVAssert(m.rowsize() == v.size());
    m.LDivEq(v); 
    return v; 
  }

  template <class T> inline const VectorView<CT>& operator/=(
      const VectorView<CT>& v, const GenMatrix<T>& m)
  {
    TMVAssert(m.IsSquare());
    TMVAssert(m.rowsize() == v.size());
    m.LDivEq(v); 
    return v; 
  }

  template <class T> inline const VectorView<T>& operator%=(
      const VectorView<T>& v, const GenMatrix<T>& m)
  {
    TMVAssert(m.IsSquare());
    TMVAssert(m.rowsize() == v.size());
    m.RDivEq(v); 
    return v; 
  }

  template <class T> inline const VectorView<CT>& operator%=(
      const VectorView<CT>& v, const GenMatrix<T>& m)
  { 
    TMVAssert(m.IsSquare());
    TMVAssert(m.rowsize() == v.size());
    m.RDivEq(v); 
    return v; 
  }

  template <class T, class Tm> inline const VectorView<T>& operator*=(
      const VectorView<T>& v, const QuotXM<T,Tm>& qxm)
  {
    TMVAssert(qxm.GetM().IsSquare());
    TMVAssert(qxm.GetM().rowsize() == v.size());
    qxm.GetM().RDivEq(v); 
    v *= qxm.GetX(); 
    return v; 
  }

  template <class T> inline const VectorView<CT>& operator*=(
      const VectorView<CT>& v, const QuotXM<T,T>& qxm)
  {
    TMVAssert(qxm.GetM().IsSquare());
    TMVAssert(qxm.GetM().rowsize() == v.size());
    qxm.GetM().RDivEq(v); 
    v *= qxm.GetX(); 
    return v; 
  }

#define GENMATRIX GenMatrix
#define GENVECTOR GenVector
#define PRODXM ProdXM
#define PRODXV ProdXV
#define QUOTVM QuotVM
#define RQUOTVM RQuotVM
#define QUOTXM QuotXM
#include "TMV_AuxQuotVM.h"
#undef GENMATRIX
#undef GENVECTOR
#undef PRODXM
#undef PRODXV
#undef QUOTVM
#undef RQUOTVM
#undef QUOTXM


  //
  // Matrix / % Matrix
  //

  template <class T, class T1, class T2> class QuotMM : 
    public MatrixComposite<T>
  {
    public:
      inline QuotMM(const T _x, const GenMatrix<T1>& _m1, 
	  const GenMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert( m1.colsize() == m2.colsize() ); }
      inline size_t colsize() const { return m2.rowsize(); }
      inline size_t rowsize() const { return m1.rowsize(); }
      inline StorageType stor() const { return BaseStorOf(m1); }
      inline T GetX() const { return x; }
      inline const GenMatrix<T1>& GetM1() const { return m1; }
      inline const GenMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToM(const MatrixView<RealType(T)>& m0) const
      {
	TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
	TMVAssert(IsReal(T()));
	m2.LDiv(m1,m0);
	MultXM(x,m0);
      }
      inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
	m2.LDiv(m1,m0);
	MultXM(x,m0);
      }
    protected:
      const T x;
      const GenMatrix<T1>& m1;
      const GenMatrix<T2>& m2;
  };

  template <class T, class T1, class T2> class TransientQuotMM : 
    public QuotMM<T,T1,T2>
  {
    public :
      inline TransientQuotMM(const T x, auto_ptr<Matrix<T1,ColMajor> > m1,
	  const GenMatrix<T2>& m2) :
	QuotMM<T,T1,T2>(x,*m1,m2), m1p(m1) {}
      inline TransientQuotMM(const TransientQuotMM<T,T1,T2>& rhs) :
	QuotMM<T,T1,T2>(rhs), m1p(rhs.m1p) {}
      inline ~TransientQuotMM() {}
      inline auto_ptr<Matrix<T1,ColMajor> > GetP() { return m1p; }

    private :
      mutable auto_ptr<Matrix<T1,ColMajor> > m1p;
  };

  template <class T, class T1, class T2> class RQuotMM : 
    public MatrixComposite<T>
  {
    public:
      inline RQuotMM(const T _x, const GenMatrix<T1>& _m1,
	  const GenMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert( m1.rowsize() == m2.rowsize() ); }
      inline size_t colsize() const { return m1.colsize(); }
      inline size_t rowsize() const { return m2.colsize(); }
      inline StorageType stor() const { return BaseStorOf(m1); }
      inline T GetX() const { return x; }
      inline const GenMatrix<T1>& GetM1() const { return m1; }
      inline const GenMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToM(const MatrixView<RealType(T)>& m0) const
      {
	TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
	TMVAssert(IsReal(T()));
	m2.RDiv(m1,m0);
	MultXM(x,m0);
      }
      inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
	m2.RDiv(m1,m0);
	MultXM(x,m0);
      }
    protected:
      const T x;
      const GenMatrix<T1>& m1;
      const GenMatrix<T2>& m2;
  };

  template <class T, class T1, class T2> class TransientRQuotMM : 
    public RQuotMM<T,T1,T2>
  {
    public :
      inline TransientRQuotMM(const T x, auto_ptr<Matrix<T1,RowMajor> > m1,
	  const GenMatrix<T2>& m2) :
	RQuotMM<T,T1,T2>(x,*m1,m2), m1p(m1) {}
      inline TransientRQuotMM(const TransientRQuotMM<T,T1,T2>& rhs) :
	RQuotMM<T,T1,T2>(rhs), m1p(rhs.m1p) {}
      inline ~TransientRQuotMM() {}
      inline auto_ptr<Matrix<T1,RowMajor> > GetP() { return m1p; }

    private :
      mutable auto_ptr<Matrix<T1,RowMajor> > m1p;
  };

  template <class T> inline const MatrixView<T>& operator/=(
      const MatrixView<T>& m1, const GenMatrix<T>& m2)
  { 
    TMVAssert(m2.colsize() == m2.rowsize());
    TMVAssert(m1.colsize() == m2.rowsize());
    m2.LDivEq(m1); 
    return m1; 
  }

  template <class T> inline const MatrixView<CT>& operator/=(
      const MatrixView<CT>& m1, const GenMatrix<T>& m2)
  { 
    TMVAssert(m2.colsize() == m2.rowsize());
    TMVAssert(m1.colsize() == m2.rowsize());
    m2.LDivEq(m1); 
    return m1; 
  }

  template <class T> inline const MatrixView<T>& operator%=(
      const MatrixView<T>& m1, const GenMatrix<T>& m2)
  { 
    TMVAssert(m2.colsize() == m2.rowsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    m2.RDivEq(m1); 
    return m1; 
  }

  template <class T> inline const MatrixView<CT>& operator%=(
      const MatrixView<CT>& m1, const GenMatrix<T>& m2)
  { 
    TMVAssert(m2.colsize() == m2.rowsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    m2.RDivEq(m1); 
    return m1; 
  }

  template <class T, class Tm> inline const MatrixView<T>& operator*=(
      const MatrixView<T>& m1, const QuotXM<T,Tm>& qxm)
  { 
    TMVAssert(qxm.GetM().IsSquare());
    TMVAssert(m1.rowsize() == qxm.GetM().rowsize());
    qxm.GetM().RDivEq(m1); 
    m1 *= qxm.GetX();
    return m1; 
  }

  template <class T> inline const MatrixView<CT>& operator*=(
      const MatrixView<CT>& m1, const QuotXM<T,T>& qxm)
  { 
    TMVAssert(qxm.GetM().IsSquare());
    TMVAssert(m1.rowsize() == qxm.GetM().rowsize());
    qxm.GetM().RDivEq(m1); 
    m1 *= qxm.GetX();
    return m1; 
  }

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 GenMatrix
#define PRODXM1 ProdXM
#define PRODXM2 ProdXM
#define QUOTXM QuotXM
#define QUOTMM QuotMM
#define RQUOTMM RQuotMM
#define TQUOTMM TransientQuotMM
#define TRQUOTMM TransientRQuotMM
#include "TMV_AuxQuotMM.h"
  // Defines things like m/m, m%m, (x*m)/m, m/(x*m), etc.
#include "TMV_AuxQuotMMa.h"
#include "TMV_AuxTQuotMMa.h"
  // Defines things like (m/m)*x, -(m/m), etc.
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTXM
#undef QUOTMM
#undef RQUOTMM
#undef TQUOTMM
#undef TRQUOTMM

} // namespace tmv

#undef CT
#undef CCT
#undef VCT

#endif