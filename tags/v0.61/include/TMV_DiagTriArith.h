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


#ifndef TMV_DiagTriArith_H
#define TMV_DiagTriArith_H

#include "TMV_DiagTriArithFunc.h"

#define CT std::complex<T>

namespace tmv {

  //
  // DiagMatrix + TriMatrix
  //

  template <class T, class T1, class T2> class SumDU : 
    public UpperTriMatrixComposite<T> 
  {
    public:
      inline SumDU(T _x1, const GenDiagMatrix<T1>& _m1, 
	  T _x2, const GenUpperTriMatrix<T2>& _m2) :
	x1(_x1),m1(_m1),x2(_x2),m2(_m2)
      { TMVAssert(m1.size() == m2.size()); }
      inline size_t size() const { return m1.size(); }
      inline StorageType stor() const { return BaseStorOf(m2); }
      inline DiagType dt() const { return NonUnitDiag; }
      inline T GetX1() const { return x1; }
      inline const GenDiagMatrix<T1>& GetM1() const { return m1; }
      inline T GetX2() const { return x2; }
      inline const GenUpperTriMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToU(const UpperTriMatrixView<RealType(T)>& m0) const
      { 
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	TMVAssert(m0.dt() == dt());
	MultXM(x2,m0=m2);
	AddVV(x1,m1.diag(),m0.diag());
      }
      inline void AssignToU(const UpperTriMatrixView<ComplexType(T)>& m0) const
      { 
	TMVAssert(m0.size() == size());
	TMVAssert(m0.dt() == dt());
	MultXM(x2,m0=m2);
	AddVV(x1,m1.diag(),m0.diag());
      }
    private:
      const T x1;
      const GenDiagMatrix<T1>& m1;
      const T x2;
      const GenUpperTriMatrix<T2>& m2;
  };

  template <class T> inline const UpperTriMatrixView<T>& operator+=(
      const UpperTriMatrixView<T>& m1, const GenDiagMatrix<T>& m2) 
  {
    TMVAssert(m1.size() == m2.size());
    TMVAssert(!m1.isunit());
    AddVV(T(1),m2.diag(),m1.diag());
    return m1; 
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator+=(
      const UpperTriMatrixView<CT>& m1, const GenDiagMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(!m1.isunit());
    AddVV(T(1),m2.diag(),m1.diag());
    return m1; 
  }

  template <class T> inline const UpperTriMatrixView<T>& operator-=(
      const UpperTriMatrixView<T>& m1, const GenDiagMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(!m1.isunit());
    AddVV(T(-1),m2.diag(),m1.diag());
    return m1; 
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator-=(
      const UpperTriMatrixView<CT>& m1, const GenDiagMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(!m1.isunit());
    AddVV(T(-1),m2.diag(),m1.diag());
    return m1; 
  }

  template <class T, class T2> inline const UpperTriMatrixView<T>& operator+=(
      const UpperTriMatrixView<T>& m, const ProdXD<T,T2>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    TMVAssert(!m.isunit());
    AddVV(pxm.GetX(),pxm.GetM().diag(),m.diag());
    return m;
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator+=(
      const UpperTriMatrixView<CT>& m, const ProdXD<T,T>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    TMVAssert(!m.isunit());
    AddVV(pxm.GetX(),pxm.GetM().diag(),m.diag());
    return m;
  }

  template <class T, class T2> inline const UpperTriMatrixView<T>& operator-=(
      const UpperTriMatrixView<T>& m, const ProdXD<T,T2>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    TMVAssert(!m.isunit());
    AddVV(-pxm.GetX(),pxm.GetM().diag(),m.diag());
    return m;
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator-=(
      const UpperTriMatrixView<CT>& m, const ProdXD<T,T>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    TMVAssert(!m.isunit());
    AddVV(-pxm.GetX(),pxm.GetM().diag(),m.diag());
    return m;
  }

  template <class T, class T1, class T2> class SumDL : 
    public LowerTriMatrixComposite<T> 
  {
    public:
      inline SumDL(T _x1, const GenDiagMatrix<T1>& _m1, 
	  T _x2, const GenLowerTriMatrix<T2>& _m2) :
	x1(_x1),m1(_m1),x2(_x2),m2(_m2)
      { TMVAssert(m1.size() == m2.size()); }
      inline size_t size() const { return m1.size(); }
      inline StorageType stor() const { return BaseStorOf(m2); }
      inline DiagType dt() const { return NonUnitDiag; }
      inline T GetX1() const { return x1; }
      inline const GenDiagMatrix<T1>& GetM1() const { return m1; }
      inline T GetX2() const { return x2; }
      inline const GenLowerTriMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToL(const LowerTriMatrixView<RealType(T)>& m0) const
      { 
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	TMVAssert(m0.dt() == dt());
	MultXM(x2,m0=m2);
	AddVV(x1,m1.diag(),m0.diag());
      }
      inline void AssignToL(const LowerTriMatrixView<ComplexType(T)>& m0) const
      { 
	TMVAssert(m0.size() == size());
	TMVAssert(m0.dt() == dt());
	MultXM(x2,m0=m2);
	AddVV(x1,m1.diag(),m0.diag());
      }
    private:
      const T x1;
      const GenDiagMatrix<T1>& m1;
      const T x2;
      const GenLowerTriMatrix<T2>& m2;
  };

  template <class T> inline const LowerTriMatrixView<T>& operator+=(
      const LowerTriMatrixView<T>& m1, const GenDiagMatrix<T>& m2) 
  {
    TMVAssert(m1.size() == m2.size());
    TMVAssert(!m1.isunit());
    AddVV(T(1),m2.diag(),m1.diag());
    return m1; 
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator+=(
      const LowerTriMatrixView<CT>& m1, const GenDiagMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(!m1.isunit());
    AddVV(T(1),m2.diag(),m1.diag());
    return m1; 
  }

  template <class T> inline const LowerTriMatrixView<T>& operator-=(
      const LowerTriMatrixView<T>& m1, const GenDiagMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(!m1.isunit());
    AddVV(T(-1),m2.diag(),m1.diag());
    return m1; 
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator-=(
      const LowerTriMatrixView<CT>& m1, const GenDiagMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(!m1.isunit());
    AddVV(T(-1),m2.diag(),m1.diag());
    return m1; 
  }

  template <class T, class T2> inline const LowerTriMatrixView<T>& operator+=(
      const LowerTriMatrixView<T>& m, const ProdXD<T,T2>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    TMVAssert(!m.isunit());
    AddVV(pxm.GetX(),pxm.GetM().diag(),m.diag());
    return m;
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator+=(
      const LowerTriMatrixView<CT>& m, const ProdXD<T,T>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    TMVAssert(!m.isunit());
    AddVV(pxm.GetX(),pxm.GetM().diag(),m.diag());
    return m;
  }

  template <class T, class T2> inline const LowerTriMatrixView<T>& operator-=(
      const LowerTriMatrixView<T>& m, const ProdXD<T,T2>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    TMVAssert(!m.isunit());
    AddVV(-pxm.GetX(),pxm.GetM().diag(),m.diag());
    return m;
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator-=(
      const LowerTriMatrixView<CT>& m, const ProdXD<T,T>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    TMVAssert(!m.isunit());
    AddVV(-pxm.GetX(),pxm.GetM().diag(),m.diag());
    return m;
  }

#define SUMMM SumDU
#define GENMATRIX1 GenDiagMatrix
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM1 ProdXD
#define PRODXM2 ProdXU
#include "TMV_AuxSumMM.h"
#include "TMV_AuxSumMMb.h"
#undef SUMMM
#undef GENMATRIX2
#undef PRODXM2
#define SUMMM SumDL
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM2 ProdXL
#include "TMV_AuxSumMM.h"
#include "TMV_AuxSumMMb.h"
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

  //
  // DiagMatrix * TriMatrix
  //

  template <class T, class T1, class T2> class ProdDU : 
    public UpperTriMatrixComposite<T>
  {
    public:
      inline ProdDU(T _x, const GenDiagMatrix<T1>& _m1,
	  const GenUpperTriMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert(m1.size() == m2.size()); }
      inline size_t size() const { return m1.size(); }
      inline StorageType stor() const { return BaseStorOf(m2); }
      inline DiagType dt() const { return NonUnitDiag; }
      inline T GetX() const { return x; }
      inline const GenDiagMatrix<T1>& GetM1() const { return m1; }
      inline const GenUpperTriMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToU(const UpperTriMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	TMVAssert(m0.dt() == dt());
	MultMM<false>(x,m1,m2,m0);
      }
      inline void AssignToU(const UpperTriMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(m0.dt() == dt());
	MultMM<false>(x,m1,m2,m0);
      }
    protected:
      T x;
      const GenDiagMatrix<T1>& m1;
      const GenUpperTriMatrix<T2>& m2;
  };

  template <class T, class T1, class T2> class ProdUD : 
    public UpperTriMatrixComposite<T>
  {
    public:
      inline ProdUD(T _x, const GenUpperTriMatrix<T1>& _m1,
	  const GenDiagMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert(m1.size() == m2.size()); }
      inline size_t size() const { return m1.size(); }
      inline StorageType stor() const { return BaseStorOf(m1); }
      inline DiagType dt() const { return NonUnitDiag; }
      inline T GetX() const { return x; }
      inline const GenUpperTriMatrix<T1>& GetM1() const { return m1; }
      inline const GenDiagMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToU(const UpperTriMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	TMVAssert(m0.dt() == dt());
	MultMM<false>(x,m1,m2,m0);
      }
      inline void AssignToU(const UpperTriMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(m0.dt() == dt());
	MultMM<false>(x,m1,m2,m0);
      }
    protected:
      T x;
      const GenUpperTriMatrix<T1>& m1;
      const GenDiagMatrix<T2>& m2;
  };

  template <class T> inline const UpperTriMatrixView<T>& operator*=(
      const UpperTriMatrixView<T>& m1, const GenDiagMatrix<T>& m2)
  { 
    TMVAssert(!m1.isunit());
    TMVAssert(m1.size() == m2.size());
    MultMM<false>(T(1),m1,m2,m1); 
    return m1; 
  }

  template <class T> inline const UpperTriMatrixView<CT>& operator*=(
      const UpperTriMatrixView<CT>& m1, const GenDiagMatrix<T>& m2)
  { 
    TMVAssert(!m1.isunit());
    TMVAssert(m1.size() == m2.size());
    MultMM<false>(T(1),m1,m2,m1); 
    return m1; 
  }

  template <class T, class T1, class T2> 
    inline const UpperTriMatrixView<T>& operator+=(
	const UpperTriMatrixView<T>& m, const ProdDU<T,T1,T2>& pmm)
    {
      TMVAssert(!m.isunit());
      TMVAssert(m.size() == pmm.size());
      MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T> inline const UpperTriMatrixView<CT>& operator+=(
	const UpperTriMatrixView<CT>& m, const ProdDU<T,T,T>& pmm)
    {
      TMVAssert(!m.isunit());
      TMVAssert(m.size() == pmm.size());
      MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T, class T1, class T2> 
    inline const UpperTriMatrixView<T>& operator-=(
	const UpperTriMatrixView<T>& m, const ProdDU<T,T1,T2>& pmm)
    { 
      TMVAssert(!m.isunit());
      TMVAssert(m.size() == pmm.size());
      MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T> inline const UpperTriMatrixView<CT>& operator-=(
	const UpperTriMatrixView<CT>& m, const ProdDU<T,T,T>& pmm)
    { 
      TMVAssert(!m.isunit());
      TMVAssert(m.size() == pmm.size());
      MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T, class T1, class T2> 
    inline const UpperTriMatrixView<T>& operator+=(
	const UpperTriMatrixView<T>& m, const ProdUD<T,T1,T2>& pmm)
    { 
      TMVAssert(!m.isunit());
      TMVAssert(m.size() == pmm.size());
      MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T> inline const UpperTriMatrixView<CT>& operator+=(
	const UpperTriMatrixView<CT>& m, const ProdUD<T,T,T>& pmm)
    { 
      TMVAssert(!m.isunit());
      TMVAssert(m.size() == pmm.size());
      MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T, class T1, class T2> 
    inline const UpperTriMatrixView<T>& operator-=(
	const UpperTriMatrixView<T>& m, const ProdUD<T,T1,T2>& pmm)
    { 
      TMVAssert(!m.isunit());
      TMVAssert(m.size() == pmm.size());
      MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T> inline const UpperTriMatrixView<CT>& operator-=(
	const UpperTriMatrixView<CT>& m, const ProdUD<T,T,T>& pmm)
    { 
      TMVAssert(!m.isunit());
      TMVAssert(m.size() == pmm.size());
      MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T, class T1, class T2> class ProdDL : 
    public LowerTriMatrixComposite<T>
  {
    public:
      inline ProdDL(T _x, const GenDiagMatrix<T1>& _m1,
	  const GenLowerTriMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert(m1.size() == m2.size()); }
      inline size_t size() const { return m1.size(); }
      inline StorageType stor() const { return BaseStorOf(m2); }
      inline DiagType dt() const { return NonUnitDiag; }
      inline T GetX() const { return x; }
      inline const GenDiagMatrix<T1>& GetM1() const { return m1; }
      inline const GenLowerTriMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToL(const LowerTriMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	TMVAssert(m0.dt() == dt());
	MultMM<false>(x,m1,m2,m0);
      }
      inline void AssignToL(const LowerTriMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(m0.dt() == dt());
	MultMM<false>(x,m1,m2,m0);
      }
    protected:
      T x;
      const GenDiagMatrix<T1>& m1;
      const GenLowerTriMatrix<T2>& m2;
  };

  template <class T, class T1, class T2> class ProdLD : 
    public LowerTriMatrixComposite<T>
  {
    public:
      inline ProdLD(T _x, const GenLowerTriMatrix<T1>& _m1,
	  const GenDiagMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert(m1.size() == m2.size()); }
      inline size_t size() const { return m1.size(); }
      inline StorageType stor() const { return BaseStorOf(m1); }
      inline DiagType dt() const { return NonUnitDiag; }
      inline T GetX() const { return x; }
      inline const GenLowerTriMatrix<T1>& GetM1() const { return m1; }
      inline const GenDiagMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToL(const LowerTriMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	TMVAssert(m0.dt() == dt());
	MultMM<false>(x,m1,m2,m0);
      }
      inline void AssignToL(const LowerTriMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(m0.dt() == dt());
	MultMM<false>(x,m1,m2,m0);
      }
    protected:
      T x;
      const GenLowerTriMatrix<T1>& m1;
      const GenDiagMatrix<T2>& m2;
  };

  template <class T> inline const LowerTriMatrixView<T>& operator*=(
      const LowerTriMatrixView<T>& m1, const GenDiagMatrix<T>& m2)
  { 
    TMVAssert(!m1.isunit());
    TMVAssert(m1.size() == m2.size());
    MultMM<false>(T(1),m1,m2,m1); 
    return m1; 
  }

  template <class T> inline const LowerTriMatrixView<CT>& operator*=(
      const LowerTriMatrixView<CT>& m1, const GenDiagMatrix<T>& m2)
  { 
    TMVAssert(!m1.isunit());
    TMVAssert(m1.size() == m2.size());
    MultMM<false>(T(1),m1,m2,m1); 
    return m1; 
  }

  template <class T, class T1, class T2> 
    inline const LowerTriMatrixView<T>& operator+=(
	const LowerTriMatrixView<T>& m, const ProdDL<T,T1,T2>& pmm)
    { 
      TMVAssert(!m.isunit());
      TMVAssert(m.size() == pmm.size());
      MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T> inline const LowerTriMatrixView<CT>& operator+=(
	const LowerTriMatrixView<CT>& m, const ProdDL<T,T,T>& pmm)
    { 
      TMVAssert(!m.isunit());
      TMVAssert(m.size() == pmm.size());
      MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T, class T1, class T2> 
    inline const LowerTriMatrixView<T>& operator-=(
	const LowerTriMatrixView<T>& m, const ProdDL<T,T1,T2>& pmm)
    { 
      TMVAssert(!m.isunit());
      TMVAssert(m.size() == pmm.size());
      MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T> inline const LowerTriMatrixView<CT>& operator-=(
	const LowerTriMatrixView<CT>& m, const ProdDL<T,T,T>& pmm)
    { 
      TMVAssert(!m.isunit());
      TMVAssert(m.size() == pmm.size());
      MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T, class T1, class T2> 
    inline const LowerTriMatrixView<T>& operator+=(
	const LowerTriMatrixView<T>& m, const ProdLD<T,T1,T2>& pmm)
    { 
      TMVAssert(!m.isunit());
      TMVAssert(m.size() == pmm.size());
      MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T> inline const LowerTriMatrixView<CT>& operator+=(
	const LowerTriMatrixView<CT>& m, const ProdLD<T,T,T>& pmm)
    { 
      TMVAssert(!m.isunit());
      TMVAssert(m.size() == pmm.size());
      MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T, class T1, class T2> 
    inline const LowerTriMatrixView<T>& operator-=(
	const LowerTriMatrixView<T>& m, const ProdLD<T,T1,T2>& pmm)
    { 
      TMVAssert(!m.isunit());
      TMVAssert(m.size() == pmm.size());
      MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T> inline const LowerTriMatrixView<CT>& operator-=(
	const LowerTriMatrixView<CT>& m, const ProdLD<T,T,T>& pmm)
    { 
      TMVAssert(!m.isunit());
      TMVAssert(m.size() == pmm.size());
      MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

#define PRODMM ProdUD
#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenDiagMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXD
#include "TMV_AuxProdMM.h"
#include "TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef PRODXM1
#define PRODMM ProdLD
#define GENMATRIX1 GenLowerTriMatrix
#define PRODXM1 ProdXL
#include "TMV_AuxProdMM.h"
#include "TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
  
#define PRODMM ProdDU
#define GENMATRIX1 GenDiagMatrix
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM1 ProdXD
#define PRODXM2 ProdXU
#include "TMV_AuxProdMM.h"
#include "TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX2
#undef PRODXM2
#define PRODMM ProdDL
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM2 ProdXL
#include "TMV_AuxProdMM.h"
#include "TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
  
  
  //
  // TriMatrix / % DiagMatrix
  //

  template <class T, class T1, class T2> class QuotUD : 
    public UpperTriMatrixComposite<T>
  {
    public:
      inline QuotUD(const T _x, const GenUpperTriMatrix<T1>& _m1,
	  const GenDiagMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert( m1.size() == m2.size() ); }
      inline size_t size() const { return m1.size(); }
      inline StorageType stor() const { return BaseStorOf(m1); }
      inline DiagType dt() const { return NonUnitDiag; }
      inline T GetX() const { return x; }
      inline const GenUpperTriMatrix<T1>& GetM1() const { return m1; }
      inline const GenDiagMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToU(const UpperTriMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit());
	MultMM<false>(x,DiagMatrix<T2>(m2.Inverse()),m1,m0);
      }
      inline void AssignToU(const UpperTriMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit());
	MultMM<false>(x,DiagMatrix<T2>(m2.Inverse()),m1,m0);
      }
    protected:
      const T x;
      const GenUpperTriMatrix<T1>& m1;
      const GenDiagMatrix<T2>& m2;
  };

  template <class T, class T1, class T2> class RQuotUD : 
    public UpperTriMatrixComposite<T>
  {
    public:
      inline RQuotUD(const T _x, const GenUpperTriMatrix<T1>& _m1,
	  const GenDiagMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert( m1.size() == m2.size() ); }
      inline size_t size() const { return m1.size(); }
      inline StorageType stor() const { return BaseStorOf(m1); }
      inline DiagType dt() const { return NonUnitDiag; }
      inline T GetX() const { return x; }
      inline const GenUpperTriMatrix<T1>& GetM1() const { return m1; }
      inline const GenDiagMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToU(const UpperTriMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit());
	MultMM<false>(x,m1,DiagMatrix<T2>(m2.Inverse()),m0);
      }
      inline void AssignToU(const UpperTriMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit());
	MultMM<false>(x,m1,DiagMatrix<T2>(m2.Inverse()),m0);
      }
    protected:
      const T x;
      const GenUpperTriMatrix<T1>& m1;
      const GenDiagMatrix<T2>& m2;
  };

  template <class T> inline const UpperTriMatrixView<T>& operator/=(
      const UpperTriMatrixView<T>& m1, const GenDiagMatrix<T>& m2)
  { MultMM(T(1),DiagMatrix<T>(m2.Inverse()),m1,m1); return m1; }

  template <class T> inline const UpperTriMatrixView<CT>& operator/=(
      const UpperTriMatrixView<CT>& m1, const GenDiagMatrix<T>& m2)
  { MultMM(T(1),DiagMatrix<T>(m2.Inverse()),m1,m1); return m1; }

  template <class T> inline const UpperTriMatrixView<T>& operator%=(
      const UpperTriMatrixView<T>& m1, const GenDiagMatrix<T>& m2)
  { MultMM(T(1),m1,DiagMatrix<T>(m2.Inverse()),m1); return m1; }

  template <class T> inline const UpperTriMatrixView<CT>& operator%=(
      const UpperTriMatrixView<CT>& m1, const GenDiagMatrix<T>& m2)
  { MultMM(T(1),m1,DiagMatrix<T>(m2.Inverse()),m1); return m1; }


  template <class T, class T1, class T2> class QuotLD : 
    public LowerTriMatrixComposite<T>
  {
    public:
      inline QuotLD(const T _x, const GenLowerTriMatrix<T1>& _m1,
	  const GenDiagMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert( m1.size() == m2.size() ); }
      inline size_t size() const { return m1.size(); }
      inline StorageType stor() const { return BaseStorOf(m1); }
      inline DiagType dt() const { return NonUnitDiag; }
      inline T GetX() const { return x; }
      inline const GenLowerTriMatrix<T1>& GetM1() const { return m1; }
      inline const GenDiagMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToL(const LowerTriMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit());
	MultMM<false>(x,DiagMatrix<T2>(m2.Inverse()),m1,m0);
      }
      inline void AssignToL(const LowerTriMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit());
	MultMM<false>(x,DiagMatrix<T2>(m2.Inverse()),m1,m0);
      }
    protected:
      const T x;
      const GenLowerTriMatrix<T1>& m1;
      const GenDiagMatrix<T2>& m2;
  };

  template <class T, class T1, class T2> class RQuotLD : 
    public LowerTriMatrixComposite<T>
  {
    public:
      inline RQuotLD(const T _x, const GenLowerTriMatrix<T1>& _m1,
	  const GenDiagMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert( m1.size() == m2.size() ); }
      inline size_t size() const { return m1.size(); }
      inline StorageType stor() const { return BaseStorOf(m1); }
      inline DiagType dt() const { return NonUnitDiag; }
      inline T GetX() const { return x; }
      inline const GenLowerTriMatrix<T1>& GetM1() const { return m1; }
      inline const GenDiagMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToL(const LowerTriMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit());
	MultMM<false>(x,m1,DiagMatrix<T2>(m2.Inverse()),m0);
      }
      inline void AssignToL(const LowerTriMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(!m0.isunit());
	MultMM<false>(x,m1,DiagMatrix<T2>(m2.Inverse()),m0);
      }
    protected:
      const T x;
      const GenLowerTriMatrix<T1>& m1;
      const GenDiagMatrix<T2>& m2;
  };

  template <class T> inline const LowerTriMatrixView<T>& operator/=(
      const LowerTriMatrixView<T>& m1, const GenDiagMatrix<T>& m2)
  { MultMM(T(1),DiagMatrix<T>(m2.Inverse()),m1,m1); return m1; }

  template <class T> inline const LowerTriMatrixView<CT>& operator/=(
      const LowerTriMatrixView<CT>& m1, const GenDiagMatrix<T>& m2)
  { MultMM(T(1),DiagMatrix<T>(m2.Inverse()),m1,m1); return m1; }

  template <class T> inline const LowerTriMatrixView<T>& operator%=(
      const LowerTriMatrixView<T>& m1, const GenDiagMatrix<T>& m2)
  { MultMM(T(1),m1,DiagMatrix<T>(m2.Inverse()),m1); return m1; }

  template <class T> inline const LowerTriMatrixView<CT>& operator%=(
      const LowerTriMatrixView<CT>& m1, const GenDiagMatrix<T>& m2)
  { MultMM(T(1),m1,DiagMatrix<T>(m2.Inverse()),m1); return m1; }


#define QUOTMM QuotUD
#define QUOTXM QuotXD
#define RQUOTMM RQuotUD
#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenDiagMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXD
#include "TMV_AuxQuotMM.h"
#include "TMV_AuxQuotMMa.h"
#undef QUOTMM
#undef RQUOTMM
#undef GENMATRIX1
#undef PRODXM1
#define QUOTMM QuotLD
#define RQUOTMM RQuotLD
#define GENMATRIX1 GenLowerTriMatrix
#define PRODXM1 ProdXL
#include "TMV_AuxQuotMM.h"
#include "TMV_AuxQuotMMa.h"
#undef QUOTMM
#undef QUOTXM
#undef RQUOTMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
    
#define GENMATRIX1 GenDiagMatrix
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM1 ProdXD
#define PRODXM2 ProdXU
#define QUOTXM QuotXU
#define TQUOTMM TransientQuotMU
#define TRQUOTMM TransientRQuotMU
#include "TMV_AuxTQuotMM.h"
#undef GENMATRIX2
#undef PRODXM2
#undef QUOTXM
#undef TQUOTMM
#undef TRQUOTMM
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM2 ProdXL
#define QUOTXM QuotXL
#define TQUOTMM TransientQuotML
#define TRQUOTMM TransientRQuotML
#include "TMV_AuxTQuotMM.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTXM
#undef TQUOTMM
#undef TRQUOTMM

} // namespace tmv

#undef CT

#endif
