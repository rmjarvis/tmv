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


#ifndef TMV_TriMatrixArith_H
#define TMV_TriMatrixArith_H

#include "TMV_TriMatrixArithFunc.h"
#include "TMV_MatrixArithFunc.h"
#include "TMV_VectorArithFunc.h"

#define CT std::complex<T>
#define CCT ConjRef<std::complex<T> >
#define VCT VarConjRef<std::complex<T> >

namespace tmv {

  template <class T, class Tv> class ProdXV;
  template <class T, class Tv> class ProdXM;

  //
  // First do everything which returns an UpperTriMatrixComposite
  // or only deals with an UpperTriMatrix, and not a LowerTriMatrix
  //
  
  template <class T, DiagType D, StorageType S, IndexStyle I, class Tx>
    inline UpperTriMatrix<T,D,S,I>& operator+=(
	UpperTriMatrix<T,D,S,I>& m, const Tx& x) 
    { m.View() += x; return m; }

  template <class T, DiagType D, StorageType S, IndexStyle I, class Tx>
    inline UpperTriMatrix<T,D,S,I>& operator-=(
	UpperTriMatrix<T,D,S,I>& m, const Tx& x) 
    { m.View() -= x; return m; }

  template <class T, DiagType D, StorageType S, IndexStyle I, class Tx>
    inline UpperTriMatrix<T,D,S,I>& operator*=(
	UpperTriMatrix<T,D,S,I>& m, const Tx& x) 
    { m.View() *= x; return m; }

  template <class T, DiagType D, StorageType S, IndexStyle I, class Tx>
    inline UpperTriMatrix<T,D,S,I>& operator/=(
	UpperTriMatrix<T,D,S,I>& m, const Tx& x) 
    { m.View() /= x; return m; }

  template <class T, DiagType D, StorageType S, IndexStyle I, class Tx>
    inline UpperTriMatrix<T,D,S,I>& operator%=(
	UpperTriMatrix<T,D,S,I>& m, const Tx& x) 
    { m.View() %= x; return m; }

  template <class T, class T2> inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m1, const GenUpperTriMatrix<T2>& m2) 
  {
    TMVAssert(m1.IsSquare());
    TMVAssert(m1.colsize() == m2.size());
    UpperTriMatrixViewOf(m1,NonUnitDiag) += m2; 
    return m1; 
  }

  template <class T, class T2> inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m1, const GenUpperTriMatrix<T2>& m2) 
  {
    TMVAssert(m1.IsSquare());
    TMVAssert(m1.colsize() == m2.size());
    UpperTriMatrixViewOf(m1,NonUnitDiag) -= m2; 
    return m1; 
  }

  //
  // Scalar * TriMatrix
  //

  template <class T, class Tm> class ProdXU : 
    public UpperTriMatrixComposite<T> 
  {
    public:
      inline ProdXU(const T _x, const GenUpperTriMatrix<Tm>& _m) :
	x(_x), m(_m) {}
      inline size_t size() const { return m.size(); }
      inline StorageType stor() const { return BaseStorOf(m); }
      inline DiagType dt() const { return NonUnitDiag; }
      inline T GetX() const { return x; }
      inline const GenUpperTriMatrix<Tm>& GetM() const { return m; }
      inline void AssignToU(const UpperTriMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == dt());
	MultXM(x,m0=m);
      }
      inline void AssignToU(const UpperTriMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == dt());
	MultXM(x,m0=m);
      }
    private:
      const T x;
      const GenUpperTriMatrix<Tm>& m;
  };

  // m*=x
  template <class T> inline const UpperTriMatrixView<T>& operator*=(
      const UpperTriMatrixView<T>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    MultXM(x,m);
    return m;
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator*=(
      const UpperTriMatrixView<CT>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    MultXM(x,m);
    return m;
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator*=(
      const UpperTriMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(!m.isunit());
    MultXM(CT(x),m);
    return m;
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator*=(
      const UpperTriMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(!m.isunit());
    MultXM(CT(x),m);
    return m;
  }

  // m/=x
  template <class T> inline const UpperTriMatrixView<T>& operator/=(
      const UpperTriMatrixView<T>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    MultXM(RealType(T)(1)/x,m);
    return m;
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator/=(
      const UpperTriMatrixView<CT>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    MultXM(T(1)/x,m);
    return m;
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator/=(
      const UpperTriMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(!m.isunit());
    MultXM(T(1)/CT(x),m);
    return m;
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator/=(
      const UpperTriMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(!m.isunit());
    MultXM(T(1)/CT(x),m);
    return m;
  }

#define GENMATRIX GenUpperTriMatrix
#define PRODXM ProdXU
#include "TMV_AuxProdXM.h"
#undef GENMATRIX
#undef PRODXM

  //
  // TriMatrix + Scalar
  //

  template <class T, class Tm> class SumUX : 
    public UpperTriMatrixComposite<T>
  {
    // x1*m + x2
    public:
      inline SumUX(T _x1, const GenUpperTriMatrix<Tm>& _m, T _x2) :
	x1(_x1), m(_m), x2(_x2) {}
      inline size_t size() const { return m.size(); }
      inline StorageType stor() const { return BaseStorOf(m); }
      inline DiagType dt() const { return NonUnitDiag; }
      inline T GetX1() const { return x1; }
      inline const GenUpperTriMatrix<Tm>& GetM() const { return m; }
      inline T GetX2() const { return x2; }
      inline void AssignToU(const UpperTriMatrixView<RealType(T)>& m0) const
      { 
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit());
	MultXM(x1,m0=m);
	m0.diag().AddToAll(REAL(x2));
      }
      inline void AssignToU(const UpperTriMatrixView<ComplexType(T)>& m0) const
      { 
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit());
	MultXM(x1,m0=m);
	m0.diag().AddToAll(ComplexType(T)(x2));
      }
    private:
      const T x1;
      const GenUpperTriMatrix<Tm>& m;
      const T x2;
  };

  // m+=x
  template <class T> inline const UpperTriMatrixView<T>& operator+=(
      const UpperTriMatrixView<T>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(x);
    return m; 
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator+=(
      const UpperTriMatrixView<CT>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(CT(x));
    return m; 
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator+=(
      const UpperTriMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(CT(x));
    return m; 
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator+=(
      const UpperTriMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(CT(x));
    return m; 
  }

  // m-=x
  template <class T> inline const UpperTriMatrixView<T>& operator-=(
      const UpperTriMatrixView<T>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(-x);
    return m; 
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator-=(
      const UpperTriMatrixView<CT>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(CT(-x));
    return m; 
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator-=(
      const UpperTriMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(-CT(x));
    return m; 
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator-=(
      const UpperTriMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(-CT(x));
    return m; 
  }

#define SUMMX SumUX
#define GENMATRIX GenUpperTriMatrix
#define PRODXM ProdXU
#include "TMV_AuxSumMX.h"
#undef SUMMX
#undef GENMATRIX
#undef PRODXM

  //
  // TriMatrix + TriMatrix
  //

  template <class T, class T1, class T2> class SumUU : 
    public UpperTriMatrixComposite<T> 
  {
    public:
      inline SumUU(T _x1, const GenUpperTriMatrix<T1>& _m1, 
	  T _x2, const GenUpperTriMatrix<T2>& _m2) :
	x1(_x1),m1(_m1),x2(_x2),m2(_m2)
      { TMVAssert(m1.size() == m2.size()); }
      inline size_t size() const { return m1.size(); }
      inline StorageType stor() const { return BaseStorOf(m1); }
      inline DiagType dt() const { return NonUnitDiag; }
      inline T GetX1() const { return x1; }
      inline const GenUpperTriMatrix<T1>& GetM1() const { return m1; }
      inline T GetX2() const { return x2; }
      inline const GenUpperTriMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToU(const UpperTriMatrixView<RealType(T)>& m0) const
      { 
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == dt());
	AddMM(x1,m1,x2,m2,m0);
      }
      inline void AssignToU(const UpperTriMatrixView<ComplexType(T)>& m0) const
      { 
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == dt());
	AddMM(x1,m1,x2,m2,m0);
      }
    private:
      const T x1;
      const GenUpperTriMatrix<T1>& m1;
      const T x2;
      const GenUpperTriMatrix<T2>& m2;
  };

  template <class T> inline const UpperTriMatrixView<T>& operator+=(
      const UpperTriMatrixView<T>& m1, const GenUpperTriMatrix<T>& m2) 
  {
    TMVAssert(m1.size() == m2.size());
    TMVAssert(!m1.isunit());
    AddMM(T(1),m2,m1); 
    return m1; 
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator+=(
      const UpperTriMatrixView<CT>& m1, const GenUpperTriMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(!m1.isunit());
    AddMM(T(1),m2,m1); 
    return m1; 
  }

  template <class T> inline const UpperTriMatrixView<T>& operator-=(
      const UpperTriMatrixView<T>& m1, const GenUpperTriMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(!m1.isunit());
    AddMM(T(-1),m2,m1); 
    return m1; 
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator-=(
      const UpperTriMatrixView<CT>& m1, const GenUpperTriMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(!m1.isunit());
    AddMM(T(-1),m2,m1); 
    return m1; 
  }

  template <class T, class T2> inline const UpperTriMatrixView<T>& operator+=(
      const UpperTriMatrixView<T>& m, const ProdXU<T,T2>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    TMVAssert(!m.isunit());
    AddMM(pxm.GetX(),pxm.GetM(),m);
    return m;
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator+=(
      const UpperTriMatrixView<CT>& m, const ProdXU<T,T>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    TMVAssert(!m.isunit());
    AddMM(pxm.GetX(),pxm.GetM(),m);
    return m;
  }

  template <class T, class T2> inline const UpperTriMatrixView<T>& operator-=(
      const UpperTriMatrixView<T>& m, const ProdXU<T,T2>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    TMVAssert(!m.isunit());
    AddMM(-pxm.GetX(),pxm.GetM(),m);
    return m;
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator-=(
      const UpperTriMatrixView<CT>& m, const ProdXU<T,T>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    TMVAssert(!m.isunit());
    AddMM(-pxm.GetX(),pxm.GetM(),m);
    return m;
  }

#define SUMMM SumUU
#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXU
#include "TMV_AuxSumMM.h"
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

  //
  // TriMatrix * TriMatrix
  //

  template <class T, class T1, class T2> class ProdUU : 
    public UpperTriMatrixComposite<T>
  {
    public:
      inline ProdUU(T _x, const GenUpperTriMatrix<T1>& _m1,
	  const GenUpperTriMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert(m1.size() == m2.size()); }
      inline size_t size() const { return m1.size(); }
      inline StorageType stor() const { return BaseStorOf(m1); }
      inline DiagType dt() const 
      { return (m1.isunit()&&m2.isunit()&&x==T(1)) ?  UnitDiag : NonUnitDiag; }
      inline T GetX() const { return x; }
      inline const GenUpperTriMatrix<T1>& GetM1() const { return m1; }
      inline const GenUpperTriMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToU(const UpperTriMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == dt());
	MultMM<false>(x,m1,m2,m0);
      }
      inline void AssignToU(const UpperTriMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == dt());
	MultMM<false>(x,m1,m2,m0);
      }
    protected:
      T x;
      const GenUpperTriMatrix<T1>& m1;
      const GenUpperTriMatrix<T2>& m2;
  };

  template <class T> inline const UpperTriMatrixView<T>& operator*=(
      const UpperTriMatrixView<T>& m1, const GenUpperTriMatrix<T>& m2)
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(!m1.isunit() || m2.isunit());
    MultMM<false>(T(1),m1,m2,m1); 
    return m1; 
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator*=(
      const UpperTriMatrixView<CT>& m1, const GenUpperTriMatrix<T>& m2)
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(!m1.isunit() || m2.isunit());
    MultMM<false>(T(1),m1,m2,m1); 
    return m1; 
  }

  template <class T, class T1, class T2> 
    inline const UpperTriMatrixView<T>& operator+=(
	const UpperTriMatrixView<T>& m, const ProdUU<T,T1,T2>& pmm)
    {
      TMVAssert(m.size() == pmm.size());
      TMVAssert(!m.isunit() || pmm.isunit());
      MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T> inline const UpperTriMatrixView<CT>& operator+=(
	const UpperTriMatrixView<CT>& m, const ProdUU<T,T,T>& pmm)
    { 
      TMVAssert(m.size() == pmm.size());
      TMVAssert(!m.isunit() || pmm.isunit());
      MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T, class T1, class T2> 
    inline const UpperTriMatrixView<T>& operator-=(
	const UpperTriMatrixView<T>& m, const ProdUU<T,T1,T2>& pmm)
    { 
      TMVAssert(m.size() == pmm.size());
      TMVAssert(!m.isunit() || pmm.isunit());
      MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T> inline const UpperTriMatrixView<CT>& operator-=(
	const UpperTriMatrixView<CT>& m, const ProdUU<T,T,T>& pmm)
    { 
      TMVAssert(m.size() == pmm.size());
      TMVAssert(!m.isunit() || pmm.isunit());
      MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }


#define PRODMM ProdUU
#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXU
#include "TMV_AuxProdMM.h"
#include "TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
  
  //
  // Scalar / TriMatrix
  //

  template <class T, class Tm> class QuotXU : 
    public UpperTriMatrixComposite<T> 
  {
    public:
      inline QuotXU(const T _x, const GenUpperTriMatrix<Tm>& _m) :
	x(_x), m(_m)  {}
      inline size_t size() const { return m.size(); }
      inline StorageType stor() const { return BaseStorOf(m); }
      inline DiagType dt() const 
      { return x==T(1) ? m.dt() : NonUnitDiag; }
      inline T GetX() const { return x; }
      inline const GenUpperTriMatrix<Tm>& GetM() const { return m; }
      inline void AssignToU(const UpperTriMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == dt());
	(m0 = m).Inverse(m0);
	MultXM(x,m0);
      }
      inline void AssignToU(const UpperTriMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == dt());
	(m0 = m).Inverse(m0);
	MultXM(x,m0);
      }
    private:
      const T x;
      const GenUpperTriMatrix<Tm>& m;
  };

#define QUOTXM QuotXU
#define GENMATRIX GenUpperTriMatrix
#define PRODXM ProdXU
#include "TMV_AuxQuotXM.h"
#undef QUOTXM
#undef GENMATRIX
#undef PRODXM

  //
  // TriMatrix / % TriMatrix
  //

  template <class T, class T1, class T2> class QuotUU : 
    public UpperTriMatrixComposite<T>
  {
    public:
      inline QuotUU(const T _x, const GenUpperTriMatrix<T1>& _m1,
	  const GenUpperTriMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert( m1.size() == m2.size() ); }
      inline size_t size() const { return m1.size(); }
      inline StorageType stor() const { return BaseStorOf(m2); }
      inline DiagType dt() const 
      { return x==T(1) && m1.dt()==m2.dt() ? m1.dt() : NonUnitDiag; }
      inline T GetX() const { return x; }
      inline const GenUpperTriMatrix<T1>& GetM1() const { return m1; }
      inline const GenUpperTriMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToU(const UpperTriMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || (m0.dt() == dt() && x==T(1)));
	m2.LDiv(m1,m0);
	MultXM(x,m0);
      }
      inline void AssignToU(const UpperTriMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || (m0.dt() == dt() && x==T(1)));
	m2.LDiv(m1,m0);
	MultXM(x,m0);
      }
    protected:
      const T x;
      const GenUpperTriMatrix<T1>& m1;
      const GenUpperTriMatrix<T2>& m2;
  };

  template <class T, class T1, class T2> class RQuotUU : 
    public UpperTriMatrixComposite<T>
  {
    public:
      inline RQuotUU(const T _x, const GenUpperTriMatrix<T1>& _m1,
	  const GenUpperTriMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert( m1.size() == m2.size() ); }
      inline size_t size() const { return m1.size(); }
      inline StorageType stor() const { return BaseStorOf(m2); }
      inline DiagType dt() const 
      { return x==T(1) && m1.dt()==m2.dt() ? m1.dt() : NonUnitDiag; }
      inline T GetX() const { return x; }
      inline const GenUpperTriMatrix<T1>& GetM1() const { return m1; }
      inline const GenUpperTriMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToU(const UpperTriMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || (m0.dt() == dt() && x==T(1)));
	m2.RDiv(m1,m0);
	MultXM(x,m0);
      }
      inline void AssignToU(const UpperTriMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || (m0.dt() == dt() && x==T(1)));
	m2.RDiv(m1,m0);
	MultXM(x,m0);
      }
    protected:
      const T x;
      const GenUpperTriMatrix<T1>& m1;
      const GenUpperTriMatrix<T2>& m2;
  };

  template <class T> inline const UpperTriMatrixView<T>& operator/=(
      const UpperTriMatrixView<T>& m1, const GenUpperTriMatrix<T>& m2)
  { m2.LDivEq(m1); return m1; }

  template <class T> inline const UpperTriMatrixView<CT>& operator/=(
      const UpperTriMatrixView<CT>& m1, const GenUpperTriMatrix<T>& m2)
  { m2.LDivEq(m1); return m1; }

  template <class T> inline const UpperTriMatrixView<T>& operator%=(
      const UpperTriMatrixView<T>& m1, const GenUpperTriMatrix<T>& m2)
  { m2.RDivEq(m1); return m1; }

  template <class T> inline const UpperTriMatrixView<CT>& operator%=(
      const UpperTriMatrixView<CT>& m1, const GenUpperTriMatrix<T>& m2)
  { m2.RDivEq(m1); return m1; }

#define QUOTMM QuotUU
#define QUOTXM QuotXU
#define RQUOTMM RQuotUU
#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXU
#include "TMV_AuxQuotMM.h"
#include "TMV_AuxQuotMMa.h"
#undef QUOTMM
#undef QUOTXM
#undef RQUOTMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
    
#define GENMATRIX GenUpperTriMatrix
#define PRODXM ProdXU
#define PRODMV ProdUV
#define PRODMTV ProdLV
#define QUOTVM QuotVU
#define QUOTXM QuotXU
#define RQUOTVM RQuotVU
#define CONSTMATRIXVIEW ConstUpperTriMatrixView
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
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM1 ProdXM
#define PRODXM2 ProdXU
#define SUMMM SumMU
#define PRODMM ProdMU
#define QUOTMM QuotMU
#define QUOTXM QuotXU
#define RQUOTMM RQuotMU
#define TQUOTMM TransientQuotMU
#define TRQUOTMM TransientRQuotMU

#include "TMV_AuxMatComposite1.h"
#include "TMV_AuxSumMMb.h"

#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef SUMMM
#undef PRODMM
#undef QUOTMM
#undef QUOTXM
#undef RQUOTMM
#undef TQUOTMM
#undef TRQUOTMM


#define GENMATRIX1 GenUpperTriMatrix
#define PRODXM1 ProdXU
#define GENMATRIX2 GenMatrix
#define PRODXM2 ProdXM
#define SUMMMa SumMU
#define SUMMM SumUM
#define PRODMM ProdUM
#define QUOTXM QuotXM
#define TQUOTMM TransientQuotMM
#define TRQUOTMM TransientRQuotMM

#include "TMV_AuxMatComposite2.h"
#include "TMV_AuxTQuotMM.h"

#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef SUMMM
#undef SUMMMa
#undef PRODMM
#undef QUOTXM
#undef TQUOTMM
#undef TRQUOTMM

  //
  // Now the things which return a LowerTriMatrixComposite
  // or only deal with a LowerTriMatrix
  //

  template <class T, DiagType D, StorageType S, IndexStyle I, class Tx>
    inline LowerTriMatrix<T,D,S,I>& operator+=(
	LowerTriMatrix<T,D,S,I>& m, const Tx& x) 
    { m.View() += x; return m; }

  template <class T, DiagType D, StorageType S, IndexStyle I, class Tx>
    inline LowerTriMatrix<T,D,S,I>& operator-=(
	LowerTriMatrix<T,D,S,I>& m, const Tx& x) 
    { m.View() -= x; return m; }

  template <class T, DiagType D, StorageType S, IndexStyle I, class Tx>
    inline LowerTriMatrix<T,D,S,I>& operator*=(
	LowerTriMatrix<T,D,S,I>& m, const Tx& x) 
    { m.View() *= x; return m; }

  template <class T, DiagType D, StorageType S, IndexStyle I, class Tx>
    inline LowerTriMatrix<T,D,S,I>& operator/=(
	LowerTriMatrix<T,D,S,I>& m, const Tx& x) 
    { m.View() /= x; return m; }

  template <class T, DiagType D, StorageType S, IndexStyle I, class Tx>
    inline LowerTriMatrix<T,D,S,I>& operator%=(
	LowerTriMatrix<T,D,S,I>& m, const Tx& x) 
    { m.View() %= x; return m; }

  template <class T, class T2> inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m1, const GenLowerTriMatrix<T2>& m2) 
  {
    TMVAssert(m1.IsSquare());
    TMVAssert(m1.colsize() == m2.size());
    LowerTriMatrixViewOf(m1,NonUnitDiag) += m2; 
    return m1; 
  }

  template <class T, class T2> inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m1, const GenLowerTriMatrix<T2>& m2) 
  { 
    TMVAssert(m1.IsSquare());
    TMVAssert(m1.colsize() == m2.size());
    LowerTriMatrixViewOf(m1,NonUnitDiag) -= m2; 
    return m1; 
  }

  //
  // Scalar * TriMatrix
  //

  template <class T, class Tm> class ProdXL : 
    public LowerTriMatrixComposite<T> 
  {
    public:
      inline ProdXL(const T _x, const GenLowerTriMatrix<Tm>& _m) :
	x(_x), m(_m) {}
      inline size_t size() const { return m.size(); }
      inline StorageType stor() const { return BaseStorOf(m); }
      inline DiagType dt() const { return NonUnitDiag; }
      inline T GetX() const { return x; }
      inline const GenLowerTriMatrix<Tm>& GetM() const { return m; }
      inline void AssignToL(const LowerTriMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == dt());
	MultXM(x,m0=m);
      }
      inline void AssignToL(const LowerTriMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == dt());
	MultXM(x,m0=m);
      }
    private:
      const T x;
      const GenLowerTriMatrix<Tm>& m;
  };

  // m*=x
  template <class T> inline const LowerTriMatrixView<T>& operator*=(
      const LowerTriMatrixView<T>& m, T x) 
  { 
    if (x != T(1)) {
      TMVAssert(!m.isunit());
      MultXM(x,m);
    }
    return m;
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator*=(
      const LowerTriMatrixView<CT>& m, T x) 
  { 
    if (x != T(1)) {
      TMVAssert(!m.isunit());
      MultXM(x,m);
    }
    return m;
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator*=(
      const LowerTriMatrixView<CT>& m, CCT x) 
  { 
    if (x != CT(1)) {
      TMVAssert(!m.isunit());
      MultXM(CT(x),m);
    }
    return m;
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator*=(
      const LowerTriMatrixView<CT>& m, VCT x) 
  { 
    if (x != CT(1)) {
      TMVAssert(!m.isunit());
      MultXM(CT(x),m);
    }
    return m;
  }

  // m/=x
  template <class T> inline const LowerTriMatrixView<T>& operator/=(
      const LowerTriMatrixView<T>& m, T x) 
  { 
    if (x != T(1)) {
      TMVAssert(!m.isunit());
      MultXM(RealType(T)(1)/x,m);
    }
    return m;
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator/=(
      const LowerTriMatrixView<CT>& m, T x) 
  { 
    if (x != T(1)) {
      TMVAssert(!m.isunit());
      MultXM(T(1)/x,m);
    }
    return m;
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator/=(
      const LowerTriMatrixView<CT>& m, CCT x) 
  { 
    if (x != CT(1)) {
      TMVAssert(!m.isunit());
      MultXM(T(1)/CT(x),m);
    }
    return m;
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator/=(
      const LowerTriMatrixView<CT>& m, VCT x) 
  { 
    if (x != CT(1)) {
      TMVAssert(!m.isunit());
      MultXM(T(1)/CT(x),m);
    }
    return m;
  }

#define GENMATRIX GenLowerTriMatrix
#define PRODXM ProdXL
#include "TMV_AuxProdXM.h"
#undef GENMATRIX
#undef PRODX


  //
  // TriMatrix + Scalar
  //

  template <class T, class Tm> class SumLX : 
    public LowerTriMatrixComposite<T>
  {
    // x1*m + x2
    public:
      inline SumLX(T _x1, const GenLowerTriMatrix<Tm>& _m, T _x2) :
	x1(_x1), m(_m), x2(_x2) {}
      inline size_t size() const { return m.size(); }
      inline StorageType stor() const { return BaseStorOf(m); }
      inline DiagType dt() const { return NonUnitDiag; }
      inline T GetX1() const { return x1; }
      inline const GenLowerTriMatrix<Tm>& GetM() const { return m; }
      inline T GetX2() const { return x2; }
      inline void AssignToL(const LowerTriMatrixView<RealType(T)>& m0) const
      { 
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit());
	MultXM(x1,m0=m);
	m0.diag().AddToAll(REAL(x2));
      }
      inline void AssignToL(const LowerTriMatrixView<ComplexType(T)>& m0) const
      { 
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit());
	MultXM(x1,m0=m);
	m0.diag().AddToAll(ComplexType(T)(x2));
      }
    private:
      const T x1;
      const GenLowerTriMatrix<Tm>& m;
      const T x2;
  };

  // m+=x
  template <class T> inline const LowerTriMatrixView<T>& operator+=(
      const LowerTriMatrixView<T>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(x);
    return m; 
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator+=(
      const LowerTriMatrixView<CT>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(CT(x));
    return m; 
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator+=(
      const LowerTriMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(CT(x));
    return m; 
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator+=(
      const LowerTriMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(CT(x));
    return m; 
  }

  // m-=x
  template <class T> inline const LowerTriMatrixView<T>& operator-=(
      const LowerTriMatrixView<T>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(-x);
    return m; 
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator-=(
      const LowerTriMatrixView<CT>& m, T x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(CT(-x));
    return m; 
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator-=(
      const LowerTriMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(-CT(x));
    return m; 
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator-=(
      const LowerTriMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(!m.isunit());
    m.diag().AddToAll(-CT(x));
    return m; 
  }

#define SUMMX SumLX
#define GENMATRIX GenLowerTriMatrix
#define PRODXM ProdXL
#include "TMV_AuxSumMX.h"
#undef SUMMX
#undef GENMATRIX
#undef PRODX


  //
  // TriMatrix + TriMatrix
  // 

  template <class T, class T1, class T2> class SumLL : 
    public LowerTriMatrixComposite<T> 
  {
    public:
      inline SumLL(T _x1, const GenLowerTriMatrix<T1>& _m1, 
	  T _x2, const GenLowerTriMatrix<T2>& _m2) :
	x1(_x1),m1(_m1),x2(_x2),m2(_m2)
      { TMVAssert(m1.size() == m2.size()); }
      inline size_t size() const { return m2.size(); }
      inline StorageType stor() const { return BaseStorOf(m1); }
      inline DiagType dt() const { return NonUnitDiag; }
      inline T GetX1() const { return x1; }
      inline const GenLowerTriMatrix<T1>& GetM1() const { return m1; }
      inline T GetX2() const { return x2; }
      inline const GenLowerTriMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToL(const LowerTriMatrixView<RealType(T)>& m0) const
      { 
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == dt());
	AddMM(x1,m1,x2,m2,m0);
      }
      inline void AssignToL(const LowerTriMatrixView<ComplexType(T)>& m0) const
      { 
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == dt());
	AddMM(x1,m1,x2,m2,m0);
      }
    private:
      T x1;
      const GenLowerTriMatrix<T1>& m1;
      T x2;
      const GenLowerTriMatrix<T2>& m2;
  };

  template <class T> inline const LowerTriMatrixView<T>& operator+=(
      const LowerTriMatrixView<T>& m1, const GenLowerTriMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(!m1.isunit());
    AddMM(T(1),m2,m1); 
    return m1; 
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator+=(
      const LowerTriMatrixView<CT>& m1, const GenLowerTriMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(!m1.isunit());
    AddMM(T(1),m2,m1); 
    return m1; 
  }

  template <class T> inline const LowerTriMatrixView<T>& operator-=(
      const LowerTriMatrixView<T>& m1, const GenLowerTriMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(!m1.isunit());
    AddMM(T(-1),m2,m1); 
    return m1; 
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator-=(
      const LowerTriMatrixView<CT>& m1, const GenLowerTriMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(!m1.isunit());
    AddMM(T(-1),m2,m1); 
    return m1; 
  }

  template <class T, class T2> inline const LowerTriMatrixView<T>& operator+=(
      const LowerTriMatrixView<T>& m, const ProdXL<T,T2>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    TMVAssert(!m.isunit());
    AddMM(pxm.GetX(),pxm.GetM(),m);
    return m;
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator+=(
      const LowerTriMatrixView<CT>& m, const ProdXL<T,T>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    TMVAssert(!m.isunit());
    AddMM(pxm.GetX(),pxm.GetM(),m);
    return m;
  }

  template <class T, class T2> inline const LowerTriMatrixView<T>& operator-=(
      const LowerTriMatrixView<T>& m, const ProdXL<T,T2>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    TMVAssert(!m.isunit());
    AddMM(-pxm.GetX(),pxm.GetM(),m);
    return m;
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator-=(
      const LowerTriMatrixView<CT>& m, const ProdXL<T,T>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    TMVAssert(!m.isunit());
    AddMM(-pxm.GetX(),pxm.GetM(),m);
    return m;
  }

#define SUMMM SumLL
#define GENMATRIX1 GenLowerTriMatrix
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM1 ProdXL
#define PRODXM2 ProdXL
#include "TMV_AuxSumMM.h"
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2


  //
  // TriMatrix * TriMatrix
  //

  template <class T, class T1, class T2> class ProdLL : 
    public LowerTriMatrixComposite<T>
  {
    public:
      inline ProdLL(T _x, const GenLowerTriMatrix<T1>& _m1,
	  const GenLowerTriMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert(m2.size() == m2.size()); }
      inline size_t size() const { return m1.size(); }
      inline StorageType stor() const { return BaseStorOf(m1); }
      inline DiagType dt() const 
      { return (m1.isunit()&&m2.isunit()&&x==T(1)) ? UnitDiag : NonUnitDiag; }
      inline T GetX() const { return x; }
      inline const GenLowerTriMatrix<T1>& GetM1() const { return m1; }
      inline const GenLowerTriMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToL(const LowerTriMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == dt());
	MultMM<false>(x,m1,m2,m0);
      }
      inline void AssignToL(const LowerTriMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == dt());
	MultMM<false>(x,m1,m2,m0);
      }
    protected:
      T x;
      const GenLowerTriMatrix<T1>& m1;
      const GenLowerTriMatrix<T2>& m2;
  };

  template <class T> inline const LowerTriMatrixView<T>& operator*=(
      const LowerTriMatrixView<T>& m1, const GenLowerTriMatrix<T>& m2)
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(!m1.isunit() || m2.isunit());
    MultMM<false>(T(1),m1,m2,m1); 
    return m1; 
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator*=(
      const LowerTriMatrixView<CT>& m1, const GenLowerTriMatrix<T>& m2)
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(!m1.isunit() || m2.isunit());
    MultMM<false>(T(1),m1,m2,m1); 
    return m1; 
  }

  template <class T, class T1, class T2> 
    inline const LowerTriMatrixView<T>& operator+=(
	const LowerTriMatrixView<T>& m, const ProdLL<T,T1,T2>& pmm)
    {
      TMVAssert(m.size() == pmm.size());
      TMVAssert(!m.isunit() || pmm.isunit());
      MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T> inline const LowerTriMatrixView<CT>& operator+=(
	const LowerTriMatrixView<CT>& m, const ProdLL<T,T,T>& pmm)
    {
      TMVAssert(m.size() == pmm.size());
      TMVAssert(!m.isunit() || pmm.isunit());
      MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T, class T1, class T2> 
    inline const LowerTriMatrixView<T>& operator-=(
	const LowerTriMatrixView<T>& m, const ProdLL<T,T1,T2>& pmm)
    { 
      TMVAssert(m.size() == pmm.size());
      TMVAssert(!m.isunit() || pmm.isunit());
      MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T> inline const LowerTriMatrixView<CT>& operator-=(
	const LowerTriMatrixView<CT>& m, const ProdLL<T,T,T>& pmm)
    {
      TMVAssert(m.size() == pmm.size());
      TMVAssert(!m.isunit() || pmm.isunit());
      MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

#define PRODMM ProdLL
#define GENMATRIX1 GenLowerTriMatrix
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM1 ProdXL
#define PRODXM2 ProdXL
#include "TMV_AuxProdMM.h"
#include "TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

  //
  // Scalar / TriMatrix
  //

  template <class T, class Tm> class QuotXL : 
    public LowerTriMatrixComposite<T> 
  {
    public:
      inline QuotXL(const T _x, const GenLowerTriMatrix<Tm>& _m) :
	x(_x), m(_m)  {}
      inline size_t size() const { return m.size(); }
      inline StorageType stor() const { return BaseStorOf(m); }
      inline DiagType dt() const { return x==T(1) ? m.dt() : NonUnitDiag; }
      inline T GetX() const { return x; }
      inline const GenLowerTriMatrix<Tm>& GetM() const { return m; }
      inline void AssignToL(const LowerTriMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == dt());
	(m0=m).Inverse(m0);
	MultXM(x,m0);
      }
      inline void AssignToL(const LowerTriMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || m0.dt() == dt());
	(m0=m).Inverse(m0);
	MultXM(x,m0);
      }
    private:
      const T x;
      const GenLowerTriMatrix<Tm>& m;
  };

#define QUOTXM QuotXL
#define GENMATRIX GenLowerTriMatrix
#define PRODXM ProdXL
#include "TMV_AuxQuotXM.h"
#undef QUOTXM
#undef GENMATRIX
#undef PRODX

  // 
  // TriMatrix / TriMatrix
  //

  template <class T, class T1, class T2> class QuotLL : 
    public LowerTriMatrixComposite<T>
  {
    public:
      inline QuotLL(const T _x, const GenLowerTriMatrix<T1>& _m1,
	  const GenLowerTriMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert( m1.size() == m2.size() ); }
      inline size_t size() const { return m1.size(); }
      inline StorageType stor() const { return BaseStorOf(m2); }
      inline DiagType dt() const 
      { return x==T(1) && m1.dt()==m2.dt() ? m1.dt() : NonUnitDiag; }
      inline T GetX() const { return x; }
      inline const GenLowerTriMatrix<T1>& GetM1() const { return m1; }
      inline const GenLowerTriMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToL(const LowerTriMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || (m0.dt() == dt() && x==T(1)));
	m2.LDiv(m1,m0);
	MultXM(x,m0);
      }
      inline void AssignToL(const LowerTriMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || (m0.dt() == dt() && x==T(1)));
	m2.LDiv(m1,m0);
	MultXM(x,m0);
      }
    protected:
      const T x;
      const GenLowerTriMatrix<T1>& m1;
      const GenLowerTriMatrix<T2>& m2;
  };

  template <class T, class T1, class T2> class RQuotLL : 
    public LowerTriMatrixComposite<T>
  {
    public:
      inline RQuotLL(const T _x, const GenLowerTriMatrix<T1>& _m1,
	  const GenLowerTriMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert( m1.size() == m2.size() ); }
      inline size_t size() const { return m1.size(); }
      inline StorageType stor() const { return BaseStorOf(m2); }
      inline DiagType dt() const 
      { return x==T(1) && m1.dt()==m2.dt() ? m1.dt() : NonUnitDiag; }
      inline T GetX() const { return x; }
      inline const GenLowerTriMatrix<T1>& GetM1() const { return m1; }
      inline const GenLowerTriMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToL(const LowerTriMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || (m0.dt() == dt() && x==T(1)));
	m2.RDiv(m1,m0);
	MultXM(x,m0);
      }
      inline void AssignToL(const LowerTriMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit() || (m0.dt() == dt() && x==T(1)));
	m2.RDiv(m1,m0);
	MultXM(x,m0);
      }
    protected:
      const T x;
      const GenLowerTriMatrix<T1>& m1;
      const GenLowerTriMatrix<T2>& m2;
  };

  template <class T> inline const LowerTriMatrixView<T>& operator/=(
      const LowerTriMatrixView<T>& m1, const GenLowerTriMatrix<T>& m2)
  { m2.LDivEq(m1); return m1; }

  template <class T> inline const LowerTriMatrixView<CT>& operator/=(
      const LowerTriMatrixView<CT>& m1, const GenLowerTriMatrix<T>& m2)
  { m2.LDivEq(m1); return m1; }

  template <class T> inline const LowerTriMatrixView<T>& operator%=(
      const LowerTriMatrixView<T>& m1, const GenLowerTriMatrix<T>& m2)
  { m2.RDivEq(m1); return m1; }

  template <class T> inline const LowerTriMatrixView<CT>& operator%=(
      const LowerTriMatrixView<CT>& m1, const GenLowerTriMatrix<T>& m2)
  { m2.RDivEq(m1); return m1; }

#define QUOTMM QuotLL
#define RQUOTMM RQuotLL
#define QUOTXM QuotXL
#define GENMATRIX1 GenLowerTriMatrix
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM1 ProdXL
#define PRODXM2 ProdXL
#include "TMV_AuxQuotMM.h"
#include "TMV_AuxQuotMMa.h"
#undef QUOTMM
#undef RQUOTMM
#undef QUOTXM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

#define GENMATRIX GenLowerTriMatrix
#define PRODXM ProdXL
#define PRODMV ProdLV
#define PRODMTV ProdUV
#define QUOTVM QuotVL
#define QUOTXM QuotXL
#define RQUOTVM RQuotVL
#define CONSTMATRIXVIEW ConstLowerTriMatrixView
#include "TMV_AuxVecComposite.h"
#undef PRODMV
#undef PRODMTV
#undef QUOTVM
#undef QUOTXM
#undef RQUOTVM
#undef CONSTMATRIXVIEW
#undef PRODXM
#undef GENMATRIX

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM1 ProdXM
#define PRODXM2 ProdXL
#define SUMMM SumML
#define PRODMM ProdML
#define QUOTMM QuotML
#define QUOTXM QuotXL
#define RQUOTMM RQuotML
#define TQUOTMM TransientQuotML
#define TRQUOTMM TransientRQuotML

#include "TMV_AuxMatComposite1.h"
#include "TMV_AuxSumMMb.h"

#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef SUMMM
#undef PRODMM
#undef QUOTMM
#undef QUOTXM
#undef RQUOTMM
#undef TQUOTMM
#undef TRQUOTMM

#define GENMATRIX1 GenLowerTriMatrix
#define PRODXM1 ProdXL
#define GENMATRIX2 GenMatrix
#define PRODXM2 ProdXM
#define SUMMMa SumML
#define SUMMM SumLM
#define PRODMM ProdLM
#define QUOTXM QuotXM
#define TQUOTMM TransientQuotMM
#define TRQUOTMM TransientRQuotMM

#include "TMV_AuxMatComposite2.h"
#include "TMV_AuxTQuotMM.h"

#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef SUMMM
#undef SUMMMa
#undef PRODMM
#undef QUOTXM
#undef TQUOTMM
#undef TRQUOTMM

  //
  // Finally, the things which mix Upper and Lower
  //
 
  template <class T, class T1, class T2> class SumUL : 
    public MatrixComposite<T> 
  {
    public:
      inline SumUL(T _x1, const GenUpperTriMatrix<T1>& _m1, 
	  T _x2, const GenLowerTriMatrix<T2>& _m2) :
	x1(_x1),m1(_m1),x2(_x2),m2(_m2)
      {  TMVAssert(m1.size() == m2.size()); }
      inline size_t colsize() const { return m1.size(); }
      inline size_t rowsize() const { return m1.size(); }
      inline StorageType stor() const { return BaseStorOf(m1); }
      inline T GetX1() const { return x1; }
      inline const GenUpperTriMatrix<T1>& GetM1() const { return m1; }
      inline T GetX2() const { return x2; }
      inline const GenLowerTriMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToM(const MatrixView<RealType(T)>& m0) const
      { 
	TMVAssert(IsReal(T()));
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	AddMM(x2,m2,x1,m1,m0);
      }
      inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
      { 
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	AddMM(x2,m2,x1,m1,m0);
      }
    private:
      T x1;
      const GenUpperTriMatrix<T1>& m1;
      T x2;
      const GenLowerTriMatrix<T2>& m2;
  };

#define SUMMM SumUL
#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXL
#include "TMV_AuxSumMM.h"
#include "TMV_AuxSumMMb.h"
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

  template <class T, class T1, class T2> class ProdUL : 
    public MatrixComposite<T>
  {
    public:
      inline ProdUL(T _x, const GenUpperTriMatrix<T1>& _m1,
	  const GenLowerTriMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert( m1.size() == m2.size()); }
      inline size_t colsize() const { return m1.size(); }
      inline size_t rowsize() const { return m1.size(); }
      inline StorageType stor() const { return BaseStorOf(m1); }
      inline T GetX() const { return x; }
      inline const GenUpperTriMatrix<T1>& GetM1() const { return m1; }
      inline const GenLowerTriMatrix<T2>& GetM2() const { return m2; }
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

    protected:
      T x;
      const GenUpperTriMatrix<T1>& m1;
      const GenLowerTriMatrix<T2>& m2;
  };

  template <class T, class T1, class T2> inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m, const ProdUL<T,T1,T2>& pmm)
  {
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
    return m; 
  }

  template <class T> inline const MatrixView<CT>& operator+=(
      const MatrixView<CT>& m, const ProdUL<T,T,T>& pmm)
  {
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
    return m; 
  }

  template <class T, class T1, class T2> inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m, const ProdUL<T,T1,T2>& pmm)
  {
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
    return m; 
  }

  template <class T> inline const MatrixView<CT>& operator-=(
      const MatrixView<CT>& m, const ProdUL<T,T,T>& pmm)
  {
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
    return m; 
  }

#define PRODMM ProdUL
#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXL
#include "TMV_AuxProdMM.h"
#include "TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
  
#define TQUOTMM TransientQuotML
#define TRQUOTMM TransientRQuotML
#define QUOTXM QuotXL
#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXL
#include "TMV_AuxTQuotMM.h"
#undef TQUOTMM
#undef TRQUOTMM
#undef QUOTXM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

  template <class T, class T1, class T2> class ProdLU : 
    public MatrixComposite<T>
  {
    public:
      inline ProdLU(T _x, const GenLowerTriMatrix<T1>& _m1,
	  const GenUpperTriMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert( m1.size() == m2.size()); }
      inline size_t colsize() const { return m1.size(); }
      inline size_t rowsize() const { return m1.size(); }
      inline StorageType stor() const { return BaseStorOf(m1); }
      inline T GetX() const { return x; }
      inline const GenLowerTriMatrix<T1>& GetM1() const { return m1; }
      inline const GenUpperTriMatrix<T2>& GetM2() const { return m2; }
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
      
    protected:
      T x;
      const GenLowerTriMatrix<T1>& m1;
      const GenUpperTriMatrix<T2>& m2;
  };

  template <class T, class T1, class T2> inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m, const ProdLU<T,T1,T2>& pmm)
  {
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
    return m; 
  }

  template <class T> inline const MatrixView<T>& operator+=(
      const MatrixView<CT>& m, const ProdLU<T,T,T>& pmm)
  {
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
    return m; 
  }

  template <class T, class T1, class T2> inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m, const ProdLU<T,T1,T2>& pmm)
  {
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
    return m; 
  }

  template <class T> inline const MatrixView<T>& operator-=(
      const MatrixView<CT>& m, const ProdLU<T,T,T>& pmm)
  {
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
    return m; 
  }

#define PRODMM ProdLU
#define GENMATRIX1 GenLowerTriMatrix
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM1 ProdXL
#define PRODXM2 ProdXU
#include "TMV_AuxProdMM.h"
#include "TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

#define TQUOTMM TransientQuotMU
#define TRQUOTMM TransientRQuotMU
#define QUOTXM QuotXU
#define GENMATRIX1 GenLowerTriMatrix
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM1 ProdXL
#define PRODXM2 ProdXU
#include "TMV_AuxTQuotMM.h"
#undef TQUOTMM
#undef TRQUOTMM
#undef QUOTXM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

} // namespace tmv

#undef CT
#undef CCT
#undef VCT

#ifdef TMV_DiagMatrixArith_H
#include "TMV_DiagTriArith.h"
#endif

#ifdef TMV_BandMatrixArith_H
#include "TMV_TriBandArith.h"
#endif

#ifdef TMV_SymMatrixArith_H
#include "TMV_TriSymArith.h"
#endif

#ifdef TMV_SymBandMatrixArith_H
#include "TMV_TriSymBandArith.h"
#endif

#endif