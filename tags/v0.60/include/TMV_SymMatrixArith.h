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


#ifndef TMV_SymMatrixArith_H
#define TMV_SymMatrixArith_H

#include "TMV_SymMatrixArithFunc.h"
#include "TMV_TriMatrixArithFunc.h"
#include "TMV_MatrixArithFunc.h"
#include "TMV_VectorArithFunc.h"

#define CT std::complex<T>
#define CCT ConjRef<std::complex<T> >
#define VCT VarConjRef<std::complex<T> >

namespace tmv {

  template <class T, class Tv> class ProdXV;
  template <class T, class Tv> class ProdXM;

  template <class T, UpLoType U, StorageType S, class Tx>
    inline SymMatrix<T,U,S>& operator+=(
	SymMatrix<T,U,S>& m, const Tx& x) 
    { m.View() += x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
    inline SymMatrix<T,U,S>& operator-=(
	SymMatrix<T,U,S>& m, const Tx& x) 
    { m.View() -= x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
    inline SymMatrix<T,U,S>& operator*=(
	SymMatrix<T,U,S>& m, const Tx& x) 
    { m.View() *= x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
    inline SymMatrix<T,U,S>& operator/=(
	SymMatrix<T,U,S>& m, const Tx& x) 
    { m.View() /= x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
    inline SymMatrix<T,U,S>& operator%=(
	SymMatrix<T,U,S>& m, const Tx& x) 
    { m.View() %= x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
    inline HermMatrix<T,U,S>& operator+=(
	HermMatrix<T,U,S>& m, const Tx& x) 
    { m.View() += x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
    inline HermMatrix<T,U,S>& operator-=(
	HermMatrix<T,U,S>& m, const Tx& x) 
    { m.View() -= x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
    inline HermMatrix<T,U,S>& operator*=(
	HermMatrix<T,U,S>& m, const Tx& x) 
    { m.View() *= x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
    inline HermMatrix<T,U,S>& operator/=(
	HermMatrix<T,U,S>& m, const Tx& x) 
    { m.View() /= x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
    inline HermMatrix<T,U,S>& operator%=(
	HermMatrix<T,U,S>& m, const Tx& x) 
    { m.View() %= x; return m; }


  //
  // Scalar * SymMatrix
  //

  template <class T, class Tm> class ProdXS : 
    public SymMatrixComposite<T> 
  {
    public:
      inline ProdXS(const T _x, const GenSymMatrix<Tm>& _m) :
	x(_x), m(_m) {}
      inline size_t size() const { return m.size(); }
      inline SymType sym() const { return m.sym(); }
      inline StorageType stor() const { return BaseStorOf(m); }
      inline T GetX() const { return x; }
      inline const GenSymMatrix<Tm>& GetM() const { return m; }
      inline void AssignToM(const MatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.colsize() == m.size());
	TMVAssert(m0.rowsize() == m.size());
	UpperTriMatrixViewOf(m0) = m.UpperTri();
	if (m.size() > 0)
	  LowerTriMatrixViewOf(m0).OffDiag() = m.LowerTri().OffDiag();
	MultXM(x,m0);
      }
      inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.colsize() == m.size());
	TMVAssert(m0.rowsize() == m.size());
	UpperTriMatrixViewOf(m0) = m.UpperTri();
	if (m.size() > 0)
	  LowerTriMatrixViewOf(m0).OffDiag() = m.LowerTri().OffDiag();
	MultXM(x,m0);
      }
      inline void AssignToS(const SymMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m.issym() || IMAG(x)==RealType(T)(0) ); 
	TMVAssert(m0.size() == m.size());
	TMVAssert(IsReal(Tm()) || m0.issym() == m.issym());
	MultXM(x,m0=m);
      }
      inline void AssignToS(const SymMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m.issym() || IMAG(x)==RealType(T)(0) ); 
	TMVAssert(m0.size() == m.size());
	TMVAssert(IsReal(Tm()) || m0.issym() == m.issym());
	MultXM(x,m0=m);
      }
    private:
      const T x;
      const GenSymMatrix<Tm>& m;
  };

  // m*=x
  template <class T> inline const SymMatrixView<T>& operator*=(
      const SymMatrixView<T>& m, T x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == RealType(T)(0));
    MultXM(x,m);
    return m;
  }

  template <class T> inline const SymMatrixView<CT>& operator*=(
      const SymMatrixView<CT>& m, T x) 
  { 
    MultXM(x,m);
    return m;
  }

  template <class T> inline const SymMatrixView<CT>& operator*=(
      const SymMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == RealType(T)(0));
    MultXM(CT(x),m);
    return m;
  }

  template <class T> inline const SymMatrixView<CT>& operator*=(
      const SymMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == RealType(T)(0));
    MultXM(CT(x),m);
    return m;
  }

  // m/=x
  template <class T> inline const SymMatrixView<T>& operator/=(
      const SymMatrixView<T>& m, T x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == RealType(T)(0));
    MultXM(RealType(T)(1)/x,m);
    return m;
  }

  template <class T> inline const SymMatrixView<CT>& operator/=(
      const SymMatrixView<CT>& m, T x) 
  { 
    MultXM(T(1)/x,m);
    return m;
  }

  template <class T> inline const SymMatrixView<CT>& operator/=(
      const SymMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == RealType(T)(0));
    MultXM(T(1)/CT(x),m);
    return m;
  }

  template <class T> inline const SymMatrixView<CT>& operator/=(
      const SymMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == RealType(T)(0));
    MultXM(T(1)/CT(x),m);
    return m;
  }

#define GENMATRIX GenSymMatrix
#define PRODXM ProdXS
#include "TMV_AuxProdXM.h"
#undef GENMATRIX
#undef PRODXM

  //
  // SymMatrix + Scalar
  //

  template <class T, class Tm> class SumSX : 
    public SymMatrixComposite<T>
  {
    public:
      inline SumSX(T _x1, const GenSymMatrix<Tm>& _m, T _x2) :
	x1(_x1), m(_m), x2(_x2) {}
      inline size_t size() const { return m.size(); }
      inline SymType sym() const { return m.sym(); }
      inline StorageType stor() const { return BaseStorOf(m); }
      inline T GetX1() const { return x1; }
      inline const GenSymMatrix<Tm>& GetM() const { return m; }
      inline T GetX2() const { return x2; }
      inline void AssignToM(const MatrixView<RealType(T)>& m0) const
      { 
	TMVAssert(IsReal(T()));
	TMVAssert(m0.colsize() == m.size());
	TMVAssert(m0.rowsize() == m.size());
	MultXM(x1,m0=m);
	m0.diag().AddToAll(REAL(x2));
      } 
      inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
      { 
	TMVAssert(m0.colsize() == m.size());
	TMVAssert(m0.rowsize() == m.size());
	MultXM(x1,m0=m);
	m0.diag().AddToAll(ComplexType(T)(x2));
      } 
      inline void AssignToS(const SymMatrixView<RealType(T)>& m0) const
      { 
	TMVAssert(IsReal(T()));
	TMVAssert(IsReal(Tm()) || m0.issym() == m.issym());
	TMVAssert(m.issym() || IMAG(x1)==RealType(T)(0));
	TMVAssert(m.issym() || IMAG(x2)==RealType(T)(0));
	TMVAssert(m0.size() == m.size());
	MultXM(x1,m0=m);
	m0.diag().AddToAll(REAL(x2));
      }
      inline void AssignToS(const SymMatrixView<ComplexType(T)>& m0) const
      { 
	TMVAssert(IsReal(Tm()) || m0.issym() == m.issym());
	TMVAssert(m.issym() || IMAG(x1)==RealType(T)(0));
	TMVAssert(m.issym() || IMAG(x2)==RealType(T)(0));
	TMVAssert(m0.size() == m.size());
	if (x1 == T(1)) m0 = m;
	else m0 = x1*m;
	m0.diag().AddToAll(ComplexType(T)(x2));
      }
    private:
      const T x1;
      const GenSymMatrix<Tm>& m;
      const T x2;
  };

  // m+=x
  template <class T> inline const SymMatrixView<T>& operator+=(
      const SymMatrixView<T>& m, T x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == RealType(T)(0));
    m.diag().AddToAll(x);
    return m; 
  }

  template <class T> inline const SymMatrixView<CT>& operator+=(
      const SymMatrixView<CT>& m, T x) 
  { 
    m.diag().AddToAll(CT(x));
    return m; 
  }

  template <class T> inline const SymMatrixView<CT>& operator+=(
      const SymMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == T(0));
    m.diag().AddToAll(CT(x));
    return m; 
  }

  template <class T> inline const SymMatrixView<CT>& operator+=(
      const SymMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == T(0));
    m.diag().AddToAll(CT(x));
    return m; 
  }

  // m-=x
  template <class T> inline const SymMatrixView<T>& operator-=(
      const SymMatrixView<T>& m, T x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == RealType(T)(0));
    m.diag().AddToAll(-x);
    return m; 
  }

  template <class T> inline const SymMatrixView<CT>& operator-=(
      const SymMatrixView<CT>& m, T x) 
  { 
    m.diag().AddToAll(CT(-x));
    return m; 
  }

  template <class T> inline const SymMatrixView<CT>& operator-=(
      const SymMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == T(0));
    m.diag().AddToAll(-CT(x));
    return m; 
  }

  template <class T> inline const SymMatrixView<CT>& operator-=(
      const SymMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == T(0));
    m.diag().AddToAll(-CT(x));
    return m; 
  }

#define SUMMX SumSX
#define GENMATRIX GenSymMatrix
#define PRODXM ProdXS
#include "TMV_AuxSumMX.h"
#undef SUMMX 
#undef GENMATRIX
#undef PRODXM

  //
  // SymMatrix + SymMatrix
  //

  template <class T, class T1, class T2> class SumSS : 
    public SymMatrixComposite<T> 
  {
    public:
      inline SumSS(T _x1, const GenSymMatrix<T1>& _m1, 
	  T _x2, const GenSymMatrix<T2>& _m2) :
	x1(_x1),m1(_m1),x2(_x2),m2(_m2)
      { TMVAssert(m1.size() == m2.size()); }
      inline size_t size() const { return m1.size(); }
      inline SymType sym() const 
      { return IsReal(T1()) ? m2.sym() : m1.sym(); }
      inline StorageType stor() const { return BaseStorOf(m1); }
      inline T GetX1() const { return x1; }
      inline const GenSymMatrix<T1>& GetM1() const { return m1; }
      inline T GetX2() const { return x2; }
      inline const GenSymMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToM(const MatrixView<RealType(T)>& m0) const
      { 
	TMVAssert(IsReal(T()));
	TMVAssert(m0.colsize() == size());
	TMVAssert(m0.rowsize() == size());
	AddMM(x1,m1,x2,m2,m0);
      } 
      inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
      { 
	TMVAssert(m0.colsize() == size());
	TMVAssert(m0.rowsize() == size());
	AddMM(x1,m1,x2,m2,m0);
      } 
      inline void AssignToS(const SymMatrixView<RealType(T)>& m0) const
      { 
	TMVAssert(IsReal(T()));
	TMVAssert(IsReal(T1()) || IsReal(T2()) || m1.sym() == m2.sym());
	TMVAssert(m0.size() == m1.size());
	TMVAssert(IsReal(T1()) || m0.issym() == m1.issym());
	TMVAssert(IsReal(T2()) || m0.issym() == m2.issym());
	TMVAssert(m0.issym() || IMAG(x1) == RealType(T)(0));
	TMVAssert(m0.issym() || IMAG(x2) == RealType(T)(0));
	AddMM(x1,m1,x2,m2,m0);
      }
      inline void AssignToS(const SymMatrixView<ComplexType(T)>& m0) const
      { 
	TMVAssert(IsReal(T1()) || IsReal(T2()) || m1.sym() == m2.sym());
	TMVAssert(m0.size() == m1.size());
	TMVAssert(IsReal(T1()) || m0.issym() == m1.issym());
	TMVAssert(IsReal(T2()) || m0.issym() == m2.issym());
	TMVAssert(m0.issym() || IMAG(x1) == RealType(T)(0));
	TMVAssert(m0.issym() || IMAG(x2) == RealType(T)(0));
	AddMM(x1,m1,x2,m2,m0);
      }
    private:
      T x1;
      const GenSymMatrix<T1>& m1;
      T x2;
      const GenSymMatrix<T2>& m2;
  };

  template <class T> inline const SymMatrixView<T>& operator+=(
      const SymMatrixView<T>& m1, const GenSymMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(IsReal(T()) || m1.sym() == m2.sym());
    AddMM(T(1),m2,m1); return m1; 
  }

  template <class T> inline const SymMatrixView<CT>& operator+=(
      const SymMatrixView<CT>& m1, const GenSymMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    AddMM(T(1),m2,m1); return m1; 
  }

  template <class T> inline const SymMatrixView<T>& operator-=(
      const SymMatrixView<T>& m1, const GenSymMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(IsReal(T()) || m1.sym() == m2.sym());
    AddMM(T(-1),m2,m1); return m1; 
  }

  template <class T> inline const SymMatrixView<CT>& operator-=(
      const SymMatrixView<CT>& m1, const GenSymMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    AddMM(T(-1),m2,m1); return m1; 
  }

  template <class T, class T2> inline const SymMatrixView<T>& operator+=(
      const SymMatrixView<T>& m, const ProdXS<T,T2>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    TMVAssert(IsReal(T2()) || m.issym() == pxm.GetM().issym());
    TMVAssert(m.issym() || IMAG(pxm.GetX()) == RealType(T)(0));
    AddMM(pxm.GetX(),pxm.GetM(),m);
    return m;
  }

  template <class T> inline const SymMatrixView<CT>& operator+=(
      const SymMatrixView<CT>& m, const ProdXS<T,T>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    AddMM(pxm.GetX(),pxm.GetM(),m);
    return m;
  }

  template <class T, class T2> inline const SymMatrixView<T>& operator-=(
      const SymMatrixView<T>& m, const ProdXS<T,T2>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    TMVAssert(IsReal(T2()) || m.issym() == pxm.GetM().issym());
    TMVAssert(m.issym() || IMAG(pxm.GetX()) == RealType(T)(0));
    AddMM(-pxm.GetX(),pxm.GetM(),m);
    return m;
  }

  template <class T> inline const SymMatrixView<CT>& operator-=(
      const SymMatrixView<CT>& m, const ProdXS<T,T>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    AddMM(-pxm.GetX(),pxm.GetM(),m);
    return m;
  }

#define SUMMM SumSS
#define GENMATRIX1 GenSymMatrix
#define GENMATRIX2 GenSymMatrix
#define PRODXM1 ProdXS
#define PRODXM2 ProdXS
#include "TMV_AuxSumMM.h"
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2


  //
  // SymMatrix * SymMatrix
  //

  template <class T, class T1, class T2> class ProdSS : 
    public MatrixComposite<T>
  {
    public:
      inline ProdSS(T _x, const GenSymMatrix<T1>& _m1,
	  const GenSymMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert(m1.size() == m2.size()); }
      inline size_t colsize() const { return m1.size(); }
      inline size_t rowsize() const { return m1.size(); }
      inline StorageType stor() const { return BaseStorOf(m1); }
      inline T GetX() const { return x; }
      inline const GenSymMatrix<T1>& GetM1() const { return m1; }
      inline const GenSymMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToM(const MatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	MultMM<false>(x,m1,m2,m0);
      }
      inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	MultMM<false>(x,m1,m2,m0);
      }
    private:
      T x;
      const GenSymMatrix<T1>& m1;
      const GenSymMatrix<T2>& m2;
  };

  template <class T, class T1, class T2> 
    inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m, const ProdSS<T,T1,T2>& pmm)
    {
      TMVAssert(m.colsize() == pmm.colsize());
      TMVAssert(m.rowsize() == pmm.rowsize());
      MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T> inline const MatrixView<CT>& operator+=(
      const MatrixView<CT>& m, const ProdSS<T,T,T>& pmm)
    {
      TMVAssert(m.colsize() == pmm.colsize());
      TMVAssert(m.rowsize() == pmm.rowsize());
      MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T, class T1, class T2> 
    inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m, const ProdSS<T,T1,T2>& pmm)
    { 
      TMVAssert(m.colsize() == pmm.colsize());
      TMVAssert(m.rowsize() == pmm.rowsize());
      MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T> inline const MatrixView<CT>& operator-=(
      const MatrixView<CT>& m, const ProdSS<T,T,T>& pmm)
    { 
      TMVAssert(m.colsize() == pmm.colsize());
      TMVAssert(m.rowsize() == pmm.rowsize());
      MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

#define PRODMM ProdSS
#define GENMATRIX1 GenSymMatrix
#define GENMATRIX2 GenSymMatrix
#define PRODXM1 ProdXS
#define PRODXM2 ProdXS
#include "TMV_AuxProdMM.h"
#include "TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2


  //
  // Vector ^ Vector
  //

  // m += (x*v^v)
  template <class T, class Tv> inline const SymMatrixView<T>& operator+=(
      const SymMatrixView<T>& m0, const OProdVV<T,Tv,Tv>& opvv)
  {
    TMVAssert(m0.size() == opvv.colsize());
    TMVAssert(m0.size() == opvv.rowsize());
    TMVAssert(opvv.GetV1().SameAs(
	  m0.isherm() ? opvv.GetV2().Conjugate() : opvv.GetV2().View()));
    Rank1Update<true>(opvv.GetX(), opvv.GetV1(), m0);
    return m0;
  }

  template <class T> inline const SymMatrixView<CT>& operator+=(
      const SymMatrixView<CT>& m0, const OProdVV<T,T,T>& opvv)
  {
    TMVAssert(m0.size() == opvv.colsize());
    TMVAssert(m0.size() == opvv.rowsize());
    TMVAssert(opvv.GetV1().SameAs(opvv.GetV2()));
    Rank1Update<true>(opvv.GetX(), opvv.GetV1(), m0);
    return m0;
  }

  // m -= (x*v^v)
  template <class T, class Tv> inline const SymMatrixView<T>& operator-=(
      const SymMatrixView<T>& m0, const OProdVV<T,Tv,Tv>& opvv)
  {
    TMVAssert(m0.size() == opvv.colsize());
    TMVAssert(m0.size() == opvv.rowsize());
    TMVAssert(opvv.GetV1().SameAs(
	  m0.isherm() ? opvv.GetV2().Conjugate() : opvv.GetV2().View()));
    Rank1Update<true>(-opvv.GetX(), opvv.GetV1(), m0);
    return m0;
  }

  template <class T> inline const SymMatrixView<CT>& operator-=(
      const SymMatrixView<CT>& m0, const OProdVV<T,T,T>& opvv)
  {
    TMVAssert(m0.size() == opvv.colsize());
    TMVAssert(m0.size() == opvv.rowsize());
    TMVAssert(opvv.GetV1().SameAs(opvv.GetV2()));
    Rank1Update<true>(-opvv.GetX(), opvv.GetV1(), m0);
    return m0;
  }

  //
  // Matrix * Matrix.Transpose()
  //

  // m += (x*m*mt)
  template <class T, class Tm> inline const SymMatrixView<T>& operator+=(
      const SymMatrixView<T>& m0, const ProdMM<T,Tm,Tm>& opmm)
  {
    TMVAssert(m0.size() == opmm.colsize());
    TMVAssert(m0.size() == opmm.rowsize());
    TMVAssert(opmm.GetM1().SameAs(
	  m0.isherm() ? opmm.GetM2().Adjoint() : opmm.GetM2().Transpose()));
    RankKUpdate<true>(opmm.GetX(), opmm.GetM1(), m0);
    return m0;
  }

  template <class T> inline const SymMatrixView<CT>& operator+=(
      const SymMatrixView<CT>& m0, const ProdMM<T,T,T>& opmm)
  {
    TMVAssert(m0.size() == opmm.colsize());
    TMVAssert(m0.size() == opmm.rowsize());
    TMVAssert(opmm.GetM1().SameAs(opmm.GetM2().Transpose()));
    RankKUpdate<true>(opmm.GetX(), opmm.GetM1(), m0);
    return m0;
  }

  template <class T, class Tm> inline const SymMatrixView<T>& operator+=(
      const SymMatrixView<T>& m0, const ProdLU<T,Tm,Tm>& opmm)
  {
    TMVAssert(m0.size() == opmm.colsize());
    TMVAssert(m0.size() == opmm.rowsize());
    TMVAssert(opmm.GetM1().SameAs(
	  m0.isherm() ? opmm.GetM2().Adjoint() : opmm.GetM2().Transpose()));
    RankKUpdate<true>(opmm.GetX(), opmm.GetM1(), m0);
    return m0;
  }

  template <class T> inline const SymMatrixView<CT>& operator+=(
      const SymMatrixView<CT>& m0, const ProdLU<T,T,T>& opmm)
  {
    TMVAssert(m0.size() == opmm.colsize());
    TMVAssert(m0.size() == opmm.rowsize());
    TMVAssert(opmm.GetM1().SameAs(opmm.GetM2().Transpose()));
    RankKUpdate<true>(opmm.GetX(), opmm.GetM1(), m0);
    return m0;
  }

  template <class T, class Tm> inline const SymMatrixView<T>& operator+=(
      const SymMatrixView<T>& m0, const ProdUL<T,Tm,Tm>& opmm)
  {
    TMVAssert(m0.size() == opmm.colsize());
    TMVAssert(m0.size() == opmm.rowsize());
    TMVAssert(opmm.GetM1().SameAs(
	  m0.isherm() ? opmm.GetM2().Adjoint() : opmm.GetM2().Transpose()));
    RankKUpdate<true>(opmm.GetX(), opmm.GetM1(), m0);
    return m0;
  }

  template <class T> inline const SymMatrixView<CT>& operator+=(
      const SymMatrixView<CT>& m0, const ProdUL<T,T,T>& opmm)
  {
    TMVAssert(m0.size() == opmm.colsize());
    TMVAssert(m0.size() == opmm.rowsize());
    TMVAssert(opmm.GetM1().SameAs(opmm.GetM2().Transpose()));
    RankKUpdate<true>(opmm.GetX(), opmm.GetM1(), m0);
    return m0;
  }

  // m -= (x*m*mt)
  template <class T, class Tm> inline const SymMatrixView<T>& operator-=(
	const SymMatrixView<T>& m0, const ProdMM<T,Tm,Tm>& opmm)
    {
      TMVAssert(m0.size() == opmm.colsize());
      TMVAssert(m0.size() == opmm.rowsize());
      TMVAssert(opmm.GetM1().SameAs(
	    m0.isherm() ? opmm.GetM2().Adjoint() : opmm.GetM2().Transpose()));
      RankKUpdate<true>(-opmm.GetX(), opmm.GetM1(), m0);
      return m0;
    }

  template <class T> inline const SymMatrixView<CT>& operator-=(
      const SymMatrixView<CT>& m0, const ProdMM<T,T,T>& opmm)
  {
    TMVAssert(m0.size() == opmm.colsize());
    TMVAssert(m0.size() == opmm.rowsize());
    TMVAssert(opmm.GetM1().SameAs(opmm.GetM2().Transpose()));
    RankKUpdate<true>(-opmm.GetX(), opmm.GetM1(), m0);
    return m0;
  }

  template <class T, class Tm> inline const SymMatrixView<T>& operator-=(
	const SymMatrixView<T>& m0, const ProdUL<T,Tm,Tm>& opmm)
    {
      TMVAssert(m0.size() == opmm.colsize());
      TMVAssert(m0.size() == opmm.rowsize());
      TMVAssert(opmm.GetM1().SameAs(
	    m0.isherm() ? opmm.GetM2().Adjoint() : opmm.GetM2().Transpose()));
      RankKUpdate<true>(-opmm.GetX(), opmm.GetM1(), m0);
      return m0;
    }

  template <class T> inline const SymMatrixView<CT>& operator-=(
      const SymMatrixView<CT>& m0, const ProdUL<T,T,T>& opmm)
  {
    TMVAssert(m0.size() == opmm.colsize());
    TMVAssert(m0.size() == opmm.rowsize());
    TMVAssert(opmm.GetM1().SameAs(opmm.GetM2().Transpose()));
    RankKUpdate<true>(-opmm.GetX(), opmm.GetM1(), m0);
    return m0;
  }

  template <class T, class Tm> inline const SymMatrixView<T>& operator-=(
	const SymMatrixView<T>& m0, const ProdLU<T,Tm,Tm>& opmm)
    {
      TMVAssert(m0.size() == opmm.colsize());
      TMVAssert(m0.size() == opmm.rowsize());
      TMVAssert(opmm.GetM1().SameAs(
	    m0.isherm() ? opmm.GetM1().Adjoint() : opmm.GetM2().Transpose()));
      RankKUpdate<true>(-opmm.GetX(), opmm.GetM1(), m0);
      return m0;
    }

  template <class T> inline const SymMatrixView<CT>& operator-=(
      const SymMatrixView<CT>& m0, const ProdLU<T,T,T>& opmm)
  {
    TMVAssert(m0.size() == opmm.colsize());
    TMVAssert(m0.size() == opmm.rowsize());
    TMVAssert(opmm.GetM1().SameAs(opmm.GetM2().Transpose()));
    RankKUpdate<true>(-opmm.GetX(), opmm.GetM1(), m0);
    return m0;
  }



#define GENMATRIX GenSymMatrix
#define GENMATRIX1 GenMatrix
#define GENMATRIX2 GenSymMatrix
#define PRODXM ProdXS
#define PRODXM1 ProdXM
#define PRODXM2 ProdXS
#define SUMMM SumMS
#define PRODMM ProdMS
#define QUOTXM QuotXS
#define QUOTMM QuotMS
#define RQUOTMM RQuotMS
#define TQUOTMM TransientQuotMS
#define TRQUOTMM TransientRQuotMS

#include "TMV_AuxMatComposite1.h"
#include "TMV_AuxSumMMb.h"
#include "TMV_AuxMatComposite3.h"

#define PRODMV ProdSV
#define QUOTVM QuotVS
#define RQUOTVM RQuotVS
#define CONSTMATRIXVIEW ConstSymMatrixView
#include "TMV_AuxVecComposite.h"
#undef PRODMV
#undef QUOTVM
#undef RQUOTVM
#undef CONSTMATRIXVIEW
#undef PRODXM
#undef GENMATRIX

  // S/S -> Matrix.  Use TransientQuot
#undef GENMATRIX1
#undef PRODXM1
#define GENMATRIX1 GenSymMatrix
#define PRODXM1 ProdXS
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
#define SUMMMa SumMS
#define SUMMM SumSM
#define PRODMM ProdSM
#define TQUOTMM TransientQuotMM
#define TRQUOTMM TransientRQuotMM
#define QUOTXM QuotXM

#include "TMV_AuxMatComposite2.h"
#include "TMV_AuxTQuotMM.h"

#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

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
#include "TMV_DiagSymArith.h"
#endif

#ifdef TMV_TriMatrixArith_H
#include "TMV_TriSymArith.h"
#endif

#ifdef TMV_BandMatrixArith_H
#include "TMV_BandSymArith.h"
#endif

#ifdef TMV_SymBandMatrixArith_H
#include "TMV_SymSymBandArith.h"
#endif

#endif
