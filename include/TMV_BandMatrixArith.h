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


#ifndef TMV_BandMatrixArith_H
#define TMV_BandMatrixArith_H

#include "TMV_BaseBandMatrix.h"
#include "TMV_BandMatrixArithFunc.h"
#include "TMV_VectorArith.h"
#include "TMV_MatrixArith.h"

#define CT std::complex<T>
#define CCT ConjRef<std::complex<T> >
#define VCT VarConjRef<std::complex<T> >

namespace tmv {

  template <class T, class Tv> class ProdXV;
  template <class T, class Tv> class ProdXM;

  template <class T, StorageType S, class Tx> 
    inline BandMatrix<T,S>& operator+=(BandMatrix<T,S>& m, const Tx& x) 
    { m.View() += x; return m; }

  template <class T, StorageType S, class Tx> 
    inline BandMatrix<T,S>& operator-=(BandMatrix<T,S>& m, const Tx& x) 
    { m.View() -= x; return m; }

  template <class T, StorageType S, class Tx> 
    inline BandMatrix<T,S>& operator*=(BandMatrix<T,S>& m, const Tx& x) 
    { m.View() *= x; return m; }

  template <class T, StorageType S, class Tx> 
    inline BandMatrix<T,S>& operator/=(BandMatrix<T,S>& m, const Tx& x) 
    { m.View() /= x; return m; }

  template <class T, StorageType S, class Tx>
    inline BandMatrix<T,S>& operator%=(BandMatrix<T,S>& m, const Tx& x) 
    { m.View() %= x; return m; }

  template <class T, class T2> inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m1, const GenBandMatrix<T2>& m2)
  { 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    BandMatrixView<T>(m1,m2.nlo(),m2.nhi()) += m2;
    return m1;
  }

  template <class T, class T2> inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m1, const GenBandMatrix<T2>& m2)
  { 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    BandMatrixView<T>(m1,m2.nlo(),m2.nhi()) -= m2;
    return m1;
  }

  //
  // Scalar * BandMatrix
  //

  template <class T, class Tm> class ProdXB : 
    public BandMatrixComposite<T> 
  {
    public:
      inline ProdXB(const T _x, const GenBandMatrix<Tm>& _m) :
	x(_x), m(_m) {}
      inline size_t colsize() const { return m.colsize(); }
      inline size_t rowsize() const { return m.rowsize(); }
      inline int nlo() const { return m.nlo(); }
      inline int nhi() const { return m.nhi(); }
      inline StorageType stor() const { return BaseStorOf(m); }
      inline T GetX() const { return x; }
      inline const GenBandMatrix<Tm>& GetM() const { return m; }
      inline void AssignToB(const BandMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	TMVAssert(m0.nlo() >= nlo());
	TMVAssert(m0.nhi() >= nhi());
	MultXM(x,m0=m);
      }
      inline void AssignToB(const BandMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	TMVAssert(m0.nlo() >= nlo());
	TMVAssert(m0.nhi() >= nhi());
	MultXM(x,m0=m);
      }
    private:
      const T x;
      const GenBandMatrix<Tm>& m;
  };

  // m*=x
  template <class T> inline const BandMatrixView<T>& operator*=(
      const BandMatrixView<T>& m, T x) 
  { MultXM(x,m); return m; }

  template <class T> inline const BandMatrixView<CT>& operator*=(
      const BandMatrixView<CT>& m, T x) 
  { MultXM(x,m); return m; }

  template <class T> inline const BandMatrixView<CT>& operator*=(
      const BandMatrixView<CT>& m, CCT x) 
  { MultXM(CT(x),m); return m; }

  template <class T> inline const BandMatrixView<CT>& operator*=(
      const BandMatrixView<CT>& m, VCT x) 
  { MultXM(CT(x),m); return m; }

  // m/=x
  template <class T> inline const BandMatrixView<T>& operator/=(
      const BandMatrixView<T>& m, T x) 
  { MultXM(RealType(T)(1)/x,m); return m; }

  template <class T> inline const BandMatrixView<CT>& operator/=(
      const BandMatrixView<CT>& m, T x) 
  { MultXM(T(1)/x,m); return m; }

  template <class T> inline const BandMatrixView<CT>& operator/=(
      const BandMatrixView<CT>& m, CCT x) 
  { MultXM(T(1)/CT(x),m); return m; }

  template <class T> inline const BandMatrixView<CT>& operator/=(
      const BandMatrixView<CT>& m, VCT x) 
  { MultXM(T(1)/CT(x),m); return m; }

#define GENMATRIX GenBandMatrix
#define PRODXM ProdXB
#include "TMV_AuxProdXM.h"
#undef GENMATRIX
#undef PRODXM


  //
  // BandMatrix + Scalar
  //

  template <class T, class Tm> class SumBX : 
    public BandMatrixComposite<T>
  {
    public:
      inline SumBX(T _x1, const GenBandMatrix<Tm>& _m, T _x2) :
	x1(_x1), m(_m), x2(_x2)
	{ TMVAssert(m.IsSquare()); }
      inline size_t colsize() const { return m.colsize(); }
      inline size_t rowsize() const { return m.rowsize(); }
      inline int nlo() const { return m.nlo(); }
      inline int nhi() const { return m.nhi(); }
      inline StorageType stor() const { return BaseStorOf(m); }
      inline T GetX1() const { return x1; }
      inline const GenBandMatrix<Tm>& GetM() const { return m; }
      inline T GetX2() const { return x2; }
      inline void AssignToB(const BandMatrixView<RealType(T)>& m0) const
      { 
	TMVAssert(IsReal(T()));
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	TMVAssert(m0.nlo() >= nlo());
	TMVAssert(m0.nhi() >= nhi());
	MultXM(x1,m0=m);
	m0.diag().AddToAll(REAL(x2));
      }
      inline void AssignToB(const BandMatrixView<ComplexType(T)>& m0) const
      { 
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	TMVAssert(m0.nlo() >= nlo());
	TMVAssert(m0.nhi() >= nhi());
	MultXM(x1,m0=m);
	m0.diag().AddToAll(ComplexType(T)(x2));
      }
    private:
      const T x1;
      const GenBandMatrix<Tm>& m;
      const T x2;
  };

  // m+=x
  template <class T> inline const BandMatrixView<T>& operator+=(
      const BandMatrixView<T>& m, T x) 
  { 
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(x); return m; 
  }

  template <class T> inline const BandMatrixView<CT>& operator+=(
      const BandMatrixView<CT>& m, T x) 
  {
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(CT(x)); return m; 
  }

  template <class T> inline const BandMatrixView<CT>& operator+=(
      const BandMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(CT(x)); return m; 
  }

  template <class T> inline const BandMatrixView<CT>& operator+=(
      const BandMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(CT(x)); return m; 
  }

  // m-=x
  template <class T> inline const BandMatrixView<T>& operator-=(
      const BandMatrixView<T>& m, T x) 
  {
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(-x); return m; 
  }

  template <class T> inline const BandMatrixView<CT>& operator-=(
      const BandMatrixView<CT>& m, T x) 
  { 
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(CT(-x)); return m; 
  }

  template <class T> inline const BandMatrixView<CT>& operator-=(
      const BandMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(-CT(x)); return m; 
  }

  template <class T> inline const BandMatrixView<CT>& operator-=(
      const BandMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(m.IsSquare());
    m.diag().AddToAll(-CT(x)); return m; 
  }

#define SUMMX SumBX
#define GENMATRIX GenBandMatrix
#define PRODXM ProdXB
#include "TMV_AuxSumMX.h"
#undef SUMMX
#undef GENMATRIX
#undef PRODXM


  //
  // BandMatrix + BandMatrix
  //
  
  template <class T, class T1, class T2> class SumBB : 
    public BandMatrixComposite<T> 
  {
    public:
      inline SumBB(T _x1, const GenBandMatrix<T1>& _m1, 
	  T _x2, const GenBandMatrix<T2>& _m2) :
	x1(_x1),m1(_m1),x2(_x2),m2(_m2)
      { 
	TMVAssert(m1.rowsize() == m2.rowsize()); 
	TMVAssert(m1.colsize() == m2.colsize()); 
      }
      inline size_t colsize() const { return m1.colsize(); }
      inline size_t rowsize() const { return m1.rowsize(); }
      inline int nlo() const { return MAX(m1.nlo(),m2.nlo()); }
      inline int nhi() const { return MAX(m1.nhi(),m2.nhi()); }
      inline StorageType stor() const 
      { return m1.stor() == m2.stor() ? BaseStorOf(m1) : DiagMajor; }
      inline T GetX1() const { return x1; }
      inline const GenBandMatrix<T1>& GetM1() const { return m1; }
      inline T GetX2() const { return x2; }
      inline const GenBandMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToB(const BandMatrixView<RealType(T)>& m0) const
      { 
	TMVAssert(IsReal(T()));
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	TMVAssert(m0.nlo() >= nlo());
	TMVAssert(m0.nhi() >= nhi());
	AddMM(x1,m1,x2,m2,m0);
      }
      inline void AssignToB(const BandMatrixView<ComplexType(T)>& m0) const
      { 
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	TMVAssert(m0.nlo() >= nlo());
	TMVAssert(m0.nhi() >= nhi());
	AddMM(x1,m1,x2,m2,m0);
      }
    private:
      T x1;
      const GenBandMatrix<T1>& m1;
      T x2;
      const GenBandMatrix<T2>& m2;
  };

  template <class T> inline const BandMatrixView<T>& operator+=(
      const BandMatrixView<T>& m1, const GenBandMatrix<T>& m2) 
  { 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    TMVAssert(m1.nlo() >= m2.nlo());
    TMVAssert(m1.nhi() >= m2.nhi());
    AddMM(T(1),m2,m1); return m1; 
  }

  template <class T> inline const BandMatrixView<CT>& operator+=(
      const BandMatrixView<CT>& m1, const GenBandMatrix<T>& m2) 
  { 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    TMVAssert(m1.nlo() >= m2.nlo());
    TMVAssert(m1.nhi() >= m2.nhi());
    AddMM(T(1),m2,m1); return m1; 
  }

  template <class T> inline const BandMatrixView<T>& operator-=(
      const BandMatrixView<T>& m1, const GenBandMatrix<T>& m2) 
  { 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    TMVAssert(m1.nlo() >= m2.nlo());
    TMVAssert(m1.nhi() >= m2.nhi());
    AddMM(T(-1),m2,m1); return m1; 
  }

  template <class T> inline const BandMatrixView<CT>& operator-=(
      const BandMatrixView<CT>& m1, const GenBandMatrix<T>& m2) 
  { 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    TMVAssert(m1.nlo() >= m2.nlo());
    TMVAssert(m1.nhi() >= m2.nhi());
    AddMM(T(-1),m2,m1); return m1; 
  }

  template <class T, class T2> inline const BandMatrixView<T>& operator+=(
      const BandMatrixView<T>& m, const ProdXB<T,T2>& pxm)
  {
    TMVAssert(m.colsize() == pxm.colsize());
    TMVAssert(m.rowsize() == pxm.rowsize());
    TMVAssert(m.nlo() >= pxm.nlo());
    TMVAssert(m.nhi() >= pxm.nhi());
    AddMM(pxm.GetX(),pxm.GetM(),m);
    return m;
  }

  template <class T> inline const BandMatrixView<CT>& operator+=(
      const BandMatrixView<CT>& m, const ProdXB<T,T>& pxm)
  {
    TMVAssert(m.colsize() == pxm.colsize());
    TMVAssert(m.rowsize() == pxm.rowsize());
    TMVAssert(m.nlo() >= pxm.nlo());
    TMVAssert(m.nhi() >= pxm.nhi());
    AddMM(pxm.GetX(),pxm.GetM(),m);
    return m;
  }

  template <class T, class T2> inline const BandMatrixView<T>& operator-=(
      const BandMatrixView<T>& m, const ProdXB<T,T2>& pxm)
  {
    TMVAssert(m.colsize() == pxm.colsize());
    TMVAssert(m.rowsize() == pxm.rowsize());
    TMVAssert(m.nlo() >= pxm.nlo());
    TMVAssert(m.nhi() >= pxm.nhi());
    AddMM(-pxm.GetX(),pxm.GetM(),m);
    return m;
  }

  template <class T> inline const BandMatrixView<CT>& operator-=(
      const BandMatrixView<CT>& m, const ProdXB<T,T>& pxm)
  {
    TMVAssert(m.colsize() == pxm.colsize());
    TMVAssert(m.rowsize() == pxm.rowsize());
    TMVAssert(m.nlo() >= pxm.nlo());
    TMVAssert(m.nhi() >= pxm.nhi());
    AddMM(-pxm.GetX(),pxm.GetM(),m);
    return m;
  }

#define SUMMM SumBB
#define GENMATRIX1 GenBandMatrix
#define GENMATRIX2 GenBandMatrix
#define PRODXM1 ProdXB
#define PRODXM2 ProdXB
#include "TMV_AuxSumMM.h"
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2


  //
  // BandMatrix * BandMatrix
  //

  template <class T, class T1, class T2> class ProdBB : 
    public BandMatrixComposite<T>
  {
    public:
      inline ProdBB(T _x, const GenBandMatrix<T1>& _m1,
	  const GenBandMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert( m1.rowsize() == m2.colsize()); }
      inline size_t colsize() const { return m1.colsize(); }
      inline size_t rowsize() const { return m2.rowsize(); }
      inline int nlo() const 
      { return MIN(int(colsize()-1),m1.nlo()+m2.nlo()); }
      inline int nhi() const 
      { return MIN(int(rowsize()-1),m1.nhi()+m2.nhi()); }
      inline StorageType stor() const 
      { 
	return m1.isrm() ? RowMajor : m2.iscm() ? ColMajor :
	  (m1.iscm() && m2.isrm()) ? ColMajor : DiagMajor; 
      }
      inline T GetX() const { return x; }
      inline const GenBandMatrix<T1>& GetM1() const { return m1; }
      inline const GenBandMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToB(const BandMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	TMVAssert(m0.nlo() >= nlo());
	TMVAssert(m0.nhi() >= nhi());
	MultMM<false>(x,m1,m2,m0);
      }
      inline void AssignToB(const BandMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	TMVAssert(m0.nlo() >= nlo());
	TMVAssert(m0.nhi() >= nhi());
	MultMM<false>(x,m1,m2,m0);
      }
    private:
      T x;
      const GenBandMatrix<T1>& m1;
      const GenBandMatrix<T2>& m2;
  };

  template <class T> inline const BandMatrixView<T>& operator*=(
      const BandMatrixView<T>& m1, const GenBandMatrix<T>& m2)
  {
    TMVAssert(m2.nlo() == 0 || m1.nlo() == int(m1.colsize())-1);
    TMVAssert(m2.nhi() == 0 || m1.nhi() == int(m1.rowsize())-1);
    MultMM<false>(T(1),m1,m2,m1); 
    return m1;
  }

  template <class T> inline const BandMatrixView<CT>& operator*=(
      const BandMatrixView<CT>& m1, const GenBandMatrix<T>& m2)
  {
    TMVAssert(m2.nlo() == 0 || m1.nlo() == int(m1.colsize())-1);
    TMVAssert(m2.nhi() == 0 || m1.nhi() == int(m1.rowsize())-1);
    MultMM<false>(T(1),m1,m2,m1);
    return m1;
  }

  template <class T, class T1, class T2>
    inline const BandMatrixView<T>& operator+=(
	const BandMatrixView<T>& m, const ProdBB<T,T1,T2>& pmm)
    {
      TMVAssert(m.colsize() == pmm.colsize());
      TMVAssert(m.rowsize() == pmm.rowsize());
      TMVAssert(m.nlo() >= pmm.nlo());
      TMVAssert(m.nhi() >= pmm.nhi());
      MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T> inline const BandMatrixView<CT>& operator+=(
      const BandMatrixView<CT>& m, const ProdBB<T,T,T>& pmm)
  {
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    TMVAssert(m.nlo() >= pmm.nlo());
    TMVAssert(m.nhi() >= pmm.nhi());
    MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
    return m; 
  }

  template <class T, class T1, class T2>
    inline const BandMatrixView<T>& operator-=(
	const BandMatrixView<T>& m, const ProdBB<T,T1,T2>& pmm)
    { 
      TMVAssert(m.colsize() == pmm.colsize());
      TMVAssert(m.rowsize() == pmm.rowsize());
      TMVAssert(m.nlo() >= pmm.nlo());
      TMVAssert(m.nhi() >= pmm.nhi());
      MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T> inline const BandMatrixView<CT>& operator-=(
	const BandMatrixView<CT>& m, const ProdBB<T,T,T>& pmm)
  { 
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    TMVAssert(m.nlo() >= pmm.nlo());
    TMVAssert(m.nhi() >= pmm.nhi());
    MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
    return m; 
  }

  template <class T, class T1, class T2, class T3>
    inline const MatrixView<T>& operator+=(
	const MatrixView<T>& m, const ProdBB<T1,T2,T3>& pmm)
    {
      TMVAssert(m.colsize() == pmm.colsize());
      TMVAssert(m.rowsize() == pmm.rowsize());
      BandMatrixView<T>(m,pmm.nlo(),pmm.nhi()) += pmm; 
      return m; 
    }

  template <class T> inline const MatrixView<CT>& operator+=(
	const MatrixView<CT>& m, const ProdBB<T,T,T>& pmm)
  {
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    BandMatrixView<CT>(m,pmm.nlo(),pmm.nhi()) += pmm; 
    return m; 
  }

  template <class T, class T1, class T2, class T3>
    inline const MatrixView<T>& operator-=(
	const MatrixView<T>& m, const ProdBB<T1,T2,T3>& pmm)
    {
      TMVAssert(m.colsize() == pmm.colsize());
      TMVAssert(m.rowsize() == pmm.rowsize());
      BandMatrixView<T>(m,pmm.nlo(),pmm.nhi()) -= pmm; 
      return m; 
    }

  template <class T> inline const MatrixView<CT>& operator-=(
	const MatrixView<CT>& m, const ProdBB<T,T,T>& pmm)
  {
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    BandMatrixView<CT>(m,pmm.nlo(),pmm.nhi()) -= pmm; 
    return m;
  }

#define PRODMM ProdBB
#define GENMATRIX1 GenBandMatrix
#define GENMATRIX2 GenBandMatrix
#define PRODXM1 ProdXB
#define PRODXM2 ProdXB
#include "TMV_AuxProdMM.h"
#include "TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX GenBandMatrix
#define PRODXM ProdXB
#define PRODMV ProdBV
#define PRODMTV ProdBV
#define QUOTVM QuotVB
#define QUOTXM QuotXB
#define RQUOTVM RQuotVB
#define CONSTMATRIXVIEW ConstBandMatrixView
#include "TMV_AuxVecComposite.h"
#undef GENMATRIX
#undef PRODXM
#undef PRODMV
#undef PRODMTV
#undef QUOTVM
#undef QUOTXM
#undef RQUOTVM
#undef CONSTMATRIXVIEW

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 GenBandMatrix
#define GENMATRIX GenBandMatrix
#define PRODXM1 ProdXM
#define PRODXM2 ProdXB
#define PRODXM ProdXB
#define SUMMM SumMB
#define PRODMM ProdMB
#define QUOTXM QuotXB
#define QUOTMM QuotMB
#define RQUOTMM RQuotMB
#define TQUOTMM TransientQuotMB
#define TRQUOTMM TransientRQuotMB

#include "TMV_AuxMatComposite1.h"
#include "TMV_AuxSumMMb.h"
#include "TMV_AuxMatComposite3.h"


  // B/B -> Matrix.  Use TransientQuot
#undef GENMATRIX1
#undef PRODXM1
#define GENMATRIX1 GenBandMatrix
#define PRODXM1 ProdXB
#include "TMV_AuxTQuotMM.h"
  

#undef GENMATRIX2
#undef PRODXM2
#undef SUMMM
#undef PRODMM
#undef QUOTMM
#undef RQUOTMM
#undef TQUOTMM
#undef TRQUOTMM
#undef QUOTXM
#define GENMATRIX2 GenMatrix
#define PRODXM2 ProdXM
#define SUMMMa SumMB
#define SUMMM SumBM
#define PRODMM ProdBM
#define TQUOTMM TransientQuotMM
#define TRQUOTMM TransientRQuotMM
#define QUOTXM QuotXM

#include "TMV_AuxMatComposite2.h"
#include "TMV_AuxTQuotMM.h"

#undef GENMATRIX1
#undef GENMATRIX2
#undef GENMATRIX
#undef PRODXM1
#undef PRODXM2
#undef PRODXM
#undef SUMMM
#undef SUMMMa
#undef PRODMM
#undef TQUOTMM
#undef TRQUOTMM
#undef QUOTXM

} // namespace tmv

#undef CT
#undef CCT
#undef VCT

#ifdef TMV_DiagMatrixArith_H
#include "TMV_DiagBandArith.h"
#endif

#ifdef TMV_TriMatrixArith_H
#include "TMV_TriBandArith.h"
#endif

#ifdef TMV_SymMatrixArith_H
#include "TMV_BandSymArith.h"
#endif

#ifdef TMV_SymBandMatrixArith_H
#include "TMV_BandSymBandArith.h"
#endif

#endif
