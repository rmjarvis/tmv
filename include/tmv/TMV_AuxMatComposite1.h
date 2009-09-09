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
// (sparse) matrix that returns a Matrix and is of the form:
// (M) op (S)
// (ie. the first operand is the normal Matrix.)
//
// Need to define the following with #define statements.
// (The given definition is for a Band Matrix.  Modify as 
// appropriate for the various other matrices.)
//
// #define GENMATRIX2 GenBandMatrix
//
// #define PRODXM ProdXB
// #define SUMMM SumMB
// #define PRODMM ProdMB
// #define QUOTMM QuotMB
// #define RQUOTMM RQuotMB
// #define TQUOTMM TransientQuotMB
// #define TRQUOTMM TransientRQuotMB


//
// Matrix + Matrix
//

template <class T, class T1, class T2> 
class SUMMM : 
public MatrixComposite<T> 
{
public:
  inline SUMMM(const T _x1, const GenMatrix<T1>& _m1, 
      const T _x2, const GENMATRIX2<T2>& _m2) :
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
  inline const GENMATRIX2<T2>& GetM2() const { return m2; }
  inline void AssignToM(const MatrixView<RealType(T)>& m0) const
  { 
    TMVAssert(IsReal(T()));
    TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
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
  const GENMATRIX2<T2>& m2;
};

template <class T> 
inline const MatrixView<T>& operator+=(
    const MatrixView<T>& m1, const GENMATRIX2<T>& m2) 
{ 
  TMVAssert(m1.colsize() == m2.colsize());
  TMVAssert(m1.rowsize() == m2.rowsize());
  AddMM(T(1),m2,m1); 
  return m1; 
}

template <class T> 
inline const MatrixView<CT>& operator+=(
    const MatrixView<CT>& m1, const GENMATRIX2<T>& m2) 
{
  TMVAssert(m1.colsize() == m2.colsize());
  TMVAssert(m1.rowsize() == m2.rowsize());
  AddMM(CT(1),m2,m1); 
  return m1; 
}

template <class T> 
inline const MatrixView<T>& operator-=(
    const MatrixView<T>& m1, const GENMATRIX2<T>& m2) 
{ 
  TMVAssert(m1.colsize() == m2.colsize());
  TMVAssert(m1.rowsize() == m2.rowsize());
  AddMM(T(-1),m2,m1); 
  return m1; 
}

template <class T> 
inline const MatrixView<CT>& operator-=(
    const MatrixView<CT>& m1, const GENMATRIX2<T>& m2) 
{ 
  TMVAssert(m1.colsize() == m2.colsize());
  TMVAssert(m1.rowsize() == m2.rowsize());
  AddMM(CT(-1),m2,m1); 
  return m1; 
}

template <class T, class T2> 
inline const MatrixView<T>& operator+=(
    const MatrixView<T>& m, const PRODXM2<T,T2>& pxm)
{
  TMVAssert(m.colsize() == pxm.colsize());
  TMVAssert(m.rowsize() == pxm.rowsize());
  AddMM(pxm.GetX(),pxm.GetM(),m);
  return m;
}

template <class T> 
inline const MatrixView<CT>& operator+=(
    const MatrixView<CT>& m, const PRODXM2<T,T>& pxm)
{
  TMVAssert(m.colsize() == pxm.colsize());
  TMVAssert(m.rowsize() == pxm.rowsize());
  AddMM(CT(pxm.GetX()),pxm.GetM(),m);
  return m;
}

template <class T, class T2> 
inline const MatrixView<T>& operator-=(
    const MatrixView<T>& m, const PRODXM2<T,T2>& pxm)
{
  TMVAssert(m.colsize() == pxm.colsize());
  TMVAssert(m.rowsize() == pxm.rowsize());
  AddMM(-pxm.GetX(),pxm.GetM(),m);
  return m;
}

template <class T> 
inline const MatrixView<CT>& operator-=(
    const MatrixView<CT>& m, const PRODXM2<T,T>& pxm)
{
  TMVAssert(m.colsize() == pxm.colsize());
  TMVAssert(m.rowsize() == pxm.rowsize());
  AddMM(CT(-pxm.GetX()),pxm.GetM(),m);
  return m;
}

#include "tmv/TMV_AuxSumMM.h"
// Defines things like m+m, m-m, (x*m)-(x*m), etc.
#include "tmv/TMV_AuxSumMMa.h"
// Defines things like -(m+m), x*(m+m), (m+m)/x, etc.

//
// Matrix * Matrix
//

template <class T, class T1, class T2> 
class PRODMM : 
public MatrixComposite<T>
{
public:
  inline PRODMM(const T _x, const GenMatrix<T1>& _m1, 
      const GENMATRIX2<T2>& _m2) :
  x(_x), m1(_m1), m2(_m2)
  { TMVAssert(m1.rowsize() == m2.colsize()) ; }
  inline size_t colsize() const { return m1.colsize(); }
  inline size_t rowsize() const { return m2.rowsize(); }
  inline StorageType stor() const { return BaseStorOf(m1); }
  inline T GetX() const { return x; }
  inline const GenMatrix<T1>& GetM1() const { return m1; }
  inline const GENMATRIX2<T2>& GetM2() const { return m2; }
  inline void AssignToM(const MatrixView<RealType(T)>& m0) const
  {
    TMVAssert(IsReal(T()));
    TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
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
  const GENMATRIX2<T2>& m2;
};

template <class T> 
inline const MatrixView<T>& operator*=(
    const MatrixView<T>& m1, const GENMATRIX2<T>& m2)
{ 
  TMVAssert(m2.colsize()==m2.rowsize());
  TMVAssert(m1.rowsize()==m2.rowsize());
  MultMM<false>(T(1),m1,m2,m1); 
  return m1; 
}

template <class T> 
inline const MatrixView<CT>& operator*=(
    const MatrixView<CT>& m1, const GENMATRIX2<T>& m2)
{ 
  TMVAssert(m2.colsize()==m2.rowsize());
  TMVAssert(m1.rowsize()==m2.rowsize());
  MultMM<false>(CT(1),m1,m2,m1); 
  return m1; 
}

template <class T, class T2> 
inline const MatrixView<T>& operator*=(
    const MatrixView<T>& m1, const PRODXM2<T,T2>& pxm)
{ 
  TMVAssert(pxm.colsize()==pxm.rowsize());
  TMVAssert(m1.rowsize()==pxm.rowsize());
  MultMM<false>(pxm.GetX(),m1,pxm.GetM(),m1); 
  return m1; 
}

template <class T> 
inline const MatrixView<CT>& operator*=(
    const MatrixView<CT>& m1, const PRODXM2<T,T>& pxm)
{ 
  TMVAssert(pxm.colsize()==pxm.rowsize());
  TMVAssert(m1.rowsize()==pxm.rowsize());
  MultMM<false>(CT(pxm.GetX()),m1,pxm.GetM(),m1); 
  return m1; 
}

template <class T, class T2, class T3> 
inline const MatrixView<T>& operator+=(
    const MatrixView<T>& m, const PRODMM<T,T2,T3>& pmm)
{ 
  TMVAssert(m.colsize() == pmm.colsize());
  TMVAssert(m.rowsize() == pmm.rowsize());
  MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
  return m; 
}

template <class T> 
inline const MatrixView<CT>& operator+=(
    const MatrixView<CT>& m, const PRODMM<T,T,T>& pmm)
{ 
  TMVAssert(m.colsize() == pmm.colsize());
  TMVAssert(m.rowsize() == pmm.rowsize());
  MultMM<true>(CT(pmm.GetX()),pmm.GetM1(),pmm.GetM2(),m); 
  return m; 
}

template <class T, class T2, class T3> 
inline const MatrixView<T>& operator-=(
    const MatrixView<T>& m, const PRODMM<T,T2,T3>& pmm)
{ 
  TMVAssert(m.colsize() == pmm.colsize());
  TMVAssert(m.rowsize() == pmm.rowsize());
  MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
  return m; 
}

template <class T> 
inline const MatrixView<CT>& operator-=(
    const MatrixView<CT>& m, const PRODMM<T,T,T>& pmm)
{ 
  TMVAssert(m.colsize() == pmm.colsize());
  TMVAssert(m.rowsize() == pmm.rowsize());
  MultMM<true>(CT(-pmm.GetX()),pmm.GetM1(),pmm.GetM2(),m); 
  return m; 
}

#include "tmv/TMV_AuxProdMM.h"
// Defines things like m*m, m*(x*m), etc.
#include "tmv/TMV_AuxProdMMa.h"
// Defines things like -(m*m), x*(m*m), etc.

//
// Matrix / % Matrix
//

template <class T, class T1, class T2> 
class QUOTMM : 
public MatrixComposite<T>
{
public:
  inline QUOTMM(const T _x, const GenMatrix<T1>& _m1, 
      const GENMATRIX2<T2>& _m2) :
  x(_x), m1(_m1), m2(_m2)
  { TMVAssert( m1.colsize() == m2.colsize() ); }
  inline size_t colsize() const { return m2.rowsize(); }
  inline size_t rowsize() const { return m1.rowsize(); }
  inline StorageType stor() const { return BaseStorOf(m1); }
  inline T GetX() const { return x; }
  inline const GenMatrix<T1>& GetM1() const { return m1; }
  inline const GENMATRIX2<T2>& GetM2() const { return m2; }
  inline void AssignToM(const MatrixView<RealType(T)>& m0) const
  {
    TMVAssert(IsReal(T()));
    TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
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
  const GENMATRIX2<T2>& m2;
};

template <class T, class T1, class T2> 
class TQUOTMM : 
public QUOTMM<T,T1,T2>
{
public :
  inline TQUOTMM(const T x, auto_ptr<Matrix<T1,ColMajor> > m1,
      const GENMATRIX2<T2>& m2) :
  QUOTMM<T,T1,T2>(x,*m1,m2), m1p(m1) {}
  inline TQUOTMM(const TQUOTMM<T,T1,T2>& rhs) :
  QUOTMM<T,T1,T2>(rhs), m1p(rhs.m1p) {}
  inline ~TQUOTMM() {}
  inline auto_ptr<Matrix<T1,ColMajor> > GetP() const { return m1p; }

private :
  mutable auto_ptr<Matrix<T1,ColMajor> > m1p;
};

template <class T, class T1, class T2> 
class RQUOTMM : 
public MatrixComposite<T>
{
public:
  inline RQUOTMM(const T _x, const GenMatrix<T1>& _m1,
      const GENMATRIX2<T2>& _m2) :
  x(_x), m1(_m1), m2(_m2)
  { TMVAssert( m1.rowsize() == m2.rowsize() ); }
  inline size_t colsize() const { return m1.colsize(); }
  inline size_t rowsize() const { return m2.colsize(); }
  inline StorageType stor() const { return BaseStorOf(m1); }
  inline T GetX() const { return x; }
  inline const GenMatrix<T1>& GetM1() const { return m1; }
  inline const GENMATRIX2<T2>& GetM2() const { return m2; }
  inline void AssignToM(const MatrixView<RealType(T)>& m0) const
  {
    TMVAssert(IsReal(T()));
    TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
    m2.RDiv(m1,m0);
    if (x != T(1)) m0 *= REAL(x);
  }
  inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
  {
    TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
    m2.RDiv(m1,m0);
    if (x != T(1)) m0 *= x;
  }
protected:
  const T x;
  const GenMatrix<T1>& m1;
  const GENMATRIX2<T2>& m2;
};

template <class T, class T1, class T2> 
class TRQUOTMM : 
public RQUOTMM<T,T1,T2>
{
public :
  inline TRQUOTMM(const T x, auto_ptr<Matrix<T1,RowMajor> > m1,
      const GENMATRIX2<T2>& m2) :
  RQUOTMM<T,T1,T2>(x,*m1,m2), m1p(m1) {}
  inline TRQUOTMM(const TRQUOTMM<T,T1,T2>& rhs) :
  RQUOTMM<T,T1,T2>(rhs), m1p(rhs.m1p) {}
  inline ~TRQUOTMM() {}
  inline auto_ptr<Matrix<T1,RowMajor> > GetP() const { return m1p; }

private :
  mutable auto_ptr<Matrix<T1,RowMajor> > m1p;
};

template <class T> 
inline const MatrixView<T>& operator/=(
    const MatrixView<T>& m1, const GENMATRIX2<T>& m2)
{ 
  TMVAssert(m2.colsize() == m2.rowsize());
  TMVAssert(m1.colsize() == m2.rowsize());
  m2.LDivEq(m1); 
  return m1; 
}

template <class T> 
inline const MatrixView<CT>& operator/=(
    const MatrixView<CT>& m1, const GENMATRIX2<T>& m2)
{ 
  TMVAssert(m2.colsize() == m2.rowsize());
  TMVAssert(m1.colsize() == m2.rowsize());
  m2.LDivEq(m1); 
  return m1; 
}

template <class T> 
inline const MatrixView<T>& operator%=(
    const MatrixView<T>& m1, const GENMATRIX2<T>& m2)
{ 
  TMVAssert(m2.colsize() == m2.rowsize());
  TMVAssert(m1.rowsize() == m2.rowsize());
  m2.RDivEq(m1); 
  return m1; 
}

template <class T> 
inline const MatrixView<CT>& operator%=(
    const MatrixView<CT>& m1, const GENMATRIX2<T>& m2)
{ 
  TMVAssert(m2.colsize() == m2.rowsize());
  TMVAssert(m1.rowsize() == m2.rowsize());
  m2.RDivEq(m1); 
  return m1; 
}


template <class T, class Tm> 
inline const MatrixView<T>& operator*=(
    const MatrixView<T>& m1, const QUOTXM<T,Tm>& qxm)
{ 
  TMVAssert(qxm.GetM().IsSquare());
  TMVAssert(m1.rowsize() == qxm.GetM().rowsize());
  qxm.GetM().RDivEq(m1); 
  m1 *= qxm.GetX();
  return m1; 
}

template <class T> 
inline const MatrixView<CT>& operator*=(
    const MatrixView<CT>& m1, const QUOTXM<T,T>& qxm)
{ 
  TMVAssert(qxm.GetM().IsSquare());
  TMVAssert(m1.rowsize() == qxm.GetM().rowsize());
  qxm.GetM().RDivEq(m1); 
  m1 *= qxm.GetX();
  return m1; 
}

#include "tmv/TMV_AuxQuotMM.h"
// Defines things like m/m, m%m, (x*m)/m, m/(x*m), etc.
#include "tmv/TMV_AuxQuotMMa.h"
#include "tmv/TMV_AuxTQuotMMa.h"
// Defines things like (m/m)*x, -(m/m), etc.
