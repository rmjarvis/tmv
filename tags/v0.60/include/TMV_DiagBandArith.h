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


#ifndef TMV_DiagBandArith_H
#define TMV_DiagBandArith_H

#define CT std::complex<T>

namespace tmv {

  //
  // DiagMatrix + BandMatrix
  //

  template <class T, class T1, class T2> class SumDB : 
    public BandMatrixComposite<T> 
  {
    public:
      inline SumDB(T _x1, const GenDiagMatrix<T1>& _m1, 
	  T _x2, const GenBandMatrix<T2>& _m2) :
	x1(_x1),m1(_m1),x2(_x2),m2(_m2)
      { 
	TMVAssert(m1.size() == m2.colsize()); 
	TMVAssert(m1.size() == m2.rowsize()); 
      }
      inline size_t colsize() const { return m1.size(); }
      inline size_t rowsize() const { return m1.size(); }
      inline int nlo() const { return m2.nlo(); }
      inline int nhi() const { return m2.nhi(); }
      inline StorageType stor() const { return BaseStorOf(m2); }
      inline T GetX1() const { return x1; }
      inline const GenDiagMatrix<T1>& GetM1() const { return m1; }
      inline T GetX2() const { return x2; }
      inline const GenBandMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToB(const BandMatrixView<RealType(T)>& m0) const
      { 
	TMVAssert(IsReal(T()));
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	TMVAssert(m0.nhi() >= nhi());
	TMVAssert(m0.nlo() >= nlo());
	AddMM(x1,BandMatrixViewOf(m1),x2,m2,m0);
      }
      inline void AssignToB(const BandMatrixView<ComplexType(T)>& m0) const
      { 
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	TMVAssert(m0.nhi() >= nhi());
	TMVAssert(m0.nlo() >= nlo());
	AddMM(x1,BandMatrixViewOf(m1),x2,m2,m0);
      }
    private:
      const T x1;
      const GenDiagMatrix<T1>& m1;
      const T x2;
      const GenBandMatrix<T2>& m2;
  };

  template <class T> inline const BandMatrixView<T>& operator+=(
      const BandMatrixView<T>& m1, const GenDiagMatrix<T>& m2) 
  {
    TMVAssert(m1.colsize() == m2.size());
    TMVAssert(m1.rowsize() == m2.size());
    AddVV(T(1),m2.diag(),m1.diag());
    return m1; 
  }

  template <class T> inline const BandMatrixView<CT>& operator+=(
      const BandMatrixView<CT>& m1, const GenDiagMatrix<T>& m2) 
  { 
    TMVAssert(m1.colsize() == m2.size());
    TMVAssert(m1.rowsize() == m2.size());
    AddVV(T(1),m2.diag(),m1.diag());
    return m1; 
  }

  template <class T> inline const BandMatrixView<T>& operator-=(
      const BandMatrixView<T>& m1, const GenDiagMatrix<T>& m2) 
  { 
    TMVAssert(m1.colsize() == m2.size());
    TMVAssert(m1.rowsize() == m2.size());
    AddVV(T(-1),m2.diag(),m1.diag());
    return m1; 
  }

  template <class T> inline const BandMatrixView<CT>& operator-=(
      const BandMatrixView<CT>& m1, const GenDiagMatrix<T>& m2) 
  { 
    TMVAssert(m1.colsize() == m2.size());
    TMVAssert(m1.rowsize() == m2.size());
    AddVV(T(-1),m2.diag(),m1.diag());
    return m1; 
  }

  template <class T, class T2> inline const BandMatrixView<T>& operator+=(
      const BandMatrixView<T>& m, const ProdXD<T,T2>& pxm)
  {
    TMVAssert(m.colsize() == pxm.size());
    TMVAssert(m.rowsize() == pxm.size());
    AddVV(pxm.GetX(),pxm.GetM().diag(),m.diag());
    return m;
  }

  template <class T> inline const BandMatrixView<CT>& operator+=(
      const BandMatrixView<CT>& m, const ProdXD<T,T>& pxm)
  {
    TMVAssert(m.colsize() == pxm.size());
    TMVAssert(m.rowsize() == pxm.size());
    AddVV(pxm.GetX(),pxm.GetM().diag(),m.diag());
    return m;
  }

  template <class T, class T2> inline const BandMatrixView<T>& operator-=(
      const BandMatrixView<T>& m, const ProdXD<T,T2>& pxm)
  {
    TMVAssert(m.colsize() == pxm.size());
    TMVAssert(m.rowsize() == pxm.size());
    AddVV(-pxm.GetX(),pxm.GetM().diag(),m.diag());
    return m;
  }

  template <class T> inline const BandMatrixView<CT>& operator-=(
      const BandMatrixView<CT>& m, const ProdXD<T,T>& pxm)
  {
    TMVAssert(m.colsize() == pxm.size());
    TMVAssert(m.rowsize() == pxm.size());
    AddVV(-pxm.GetX(),pxm.GetM().diag(),m.diag());
    return m;
  }

  template <class T> inline const DiagMatrixView<T>& operator+=(
      const DiagMatrixView<T>& m1, const GenBandMatrix<T>& m2) 
  {
    TMVAssert(m2.colsize() == m1.size());
    TMVAssert(m2.rowsize() == m1.size());
    TMVAssert(m2.nlo() == 0 && m2.nhi() == 0);
    AddVV(T(1),m2.diag(),m1.diag());
    return m1; 
  }

  template <class T> inline const DiagMatrixView<CT>& operator+=(
      const DiagMatrixView<CT>& m1, const GenBandMatrix<T>& m2) 
  { 
    TMVAssert(m2.colsize() == m1.size());
    TMVAssert(m2.rowsize() == m1.size());
    TMVAssert(m2.nlo() == 0 && m2.nhi() == 0);
    AddVV(T(1),m2.diag(),m1.diag());
    return m1; 
  }

  template <class T> inline const DiagMatrixView<T>& operator-=(
      const DiagMatrixView<T>& m1, const GenBandMatrix<T>& m2) 
  { 
    TMVAssert(m2.colsize() == m1.size());
    TMVAssert(m2.rowsize() == m1.size());
    TMVAssert(m2.nlo() == 0 && m2.nhi() == 0);
    AddVV(T(-1),m2.diag(),m1.diag());
    return m1; 
  }

  template <class T> inline const DiagMatrixView<CT>& operator-=(
      const DiagMatrixView<CT>& m1, const GenBandMatrix<T>& m2) 
  { 
    TMVAssert(m2.colsize() == m1.size());
    TMVAssert(m2.rowsize() == m1.size());
    TMVAssert(m2.nlo() == 0 && m2.nhi() == 0);
    AddVV(T(-1),m2.diag(),m1.diag());
    return m1; 
  }

  template <class T, class T2> inline const DiagMatrixView<T>& operator+=(
      const DiagMatrixView<T>& m, const ProdXB<T,T2>& pxm)
  {
    TMVAssert(pxm.colsize() == m.size());
    TMVAssert(pxm.rowsize() == m.size());
    TMVAssert(pxm.nlo() == 0 && pxm.nhi() == 0);
    AddVV(pxm.GetX(),pxm.GetM().diag(),m.diag());
    return m;
  }

  template <class T> inline const DiagMatrixView<CT>& operator+=(
      const DiagMatrixView<CT>& m, const ProdXB<T,T>& pxm)
  {
    TMVAssert(pxm.colsize() == m.size());
    TMVAssert(pxm.rowsize() == m.size());
    TMVAssert(pxm.nlo() == 0 && pxm.nhi() == 0);
    AddVV(pxm.GetX(),pxm.GetM().diag(),m.diag());
    return m;
  }

  template <class T, class T2> inline const DiagMatrixView<T>& operator-=(
      const DiagMatrixView<T>& m, const ProdXB<T,T2>& pxm)
  {
    TMVAssert(pxm.colsize() == m.size());
    TMVAssert(pxm.rowsize() == m.size());
    TMVAssert(pxm.nlo() == 0 && pxm.nhi() == 0);
    AddVV(-pxm.GetX(),pxm.GetM().diag(),m.diag());
    return m;
  }

  template <class T> inline const DiagMatrixView<CT>& operator-=(
      const DiagMatrixView<CT>& m, const ProdXB<T,T>& pxm)
  {
    TMVAssert(pxm.colsize() == m.size());
    TMVAssert(pxm.rowsize() == m.size());
    TMVAssert(pxm.nlo() == 0 && pxm.nhi() == 0);
    AddVV(-pxm.GetX(),pxm.GetM().diag(),m.diag());
    return m;
  }

#define SUMMM SumDB
#define GENMATRIX1 GenDiagMatrix
#define GENMATRIX2 GenBandMatrix
#define PRODXM1 ProdXD
#define PRODXM2 ProdXB
#include "TMV_AuxSumMM.h"
#include "TMV_AuxSumMMb.h"
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

  //
  // DiagMatrix * BandMatrix
  //

  template <class T, class T1, class T2> class ProdDB : 
    public BandMatrixComposite<T>
  {
    public:
      inline ProdDB(T _x, const GenDiagMatrix<T1>& _m1,
	  const GenBandMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert(m1.size() == m2.colsize()); }
      inline size_t colsize() const { return m2.colsize(); }
      inline size_t rowsize() const { return m2.rowsize(); }
      inline int nlo() const { return m2.nlo(); }
      inline int nhi() const { return m2.nhi(); }
      inline StorageType stor() const { return BaseStorOf(m2); }
      inline T GetX() const { return x; }
      inline const GenDiagMatrix<T1>& GetM1() const { return m1; }
      inline const GenBandMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToB(const BandMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	TMVAssert(m0.nhi() >= nhi());
	TMVAssert(m0.nlo() >= nlo());
	MultMM<false>(x,BandMatrixViewOf(m1),m2,m0);
      }
      inline void AssignToB(const BandMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	TMVAssert(m0.nhi() >= nhi());
	TMVAssert(m0.nlo() >= nlo());
	MultMM<false>(x,BandMatrixViewOf(m1),m2,m0);
      }
    protected:
      T x;
      const GenDiagMatrix<T1>& m1;
      const GenBandMatrix<T2>& m2;
  };

  template <class T, class T1, class T2> class ProdBD : 
    public BandMatrixComposite<T>
  {
    public:
      inline ProdBD(T _x, const GenBandMatrix<T1>& _m1,
	  const GenDiagMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert(m1.rowsize() == m2.size()); }
      inline size_t colsize() const { return m1.colsize(); }
      inline size_t rowsize() const { return m1.rowsize(); }
      inline int nlo() const { return m1.nlo(); }
      inline int nhi() const { return m1.nhi(); }
      inline StorageType stor() const { return BaseStorOf(m1); }
      inline T GetX() const { return x; }
      inline const GenBandMatrix<T1>& GetM1() const { return m1; }
      inline const GenDiagMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToB(const BandMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	TMVAssert(m0.nhi() >= nhi());
	TMVAssert(m0.nlo() >= nlo());
	MultMM<false>(x,m1,BandMatrixViewOf(m2),m0);
      }
      inline void AssignToB(const BandMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	TMVAssert(m0.nhi() >= nhi());
	TMVAssert(m0.nlo() >= nlo());
	MultMM<false>(x,m1,BandMatrixViewOf(m2),m0);
      }
    protected:
      T x;
      const GenBandMatrix<T1>& m1;
      const GenDiagMatrix<T2>& m2;
  };

  template <class T> inline const BandMatrixView<T>& operator*=(
      const BandMatrixView<T>& m1, const GenDiagMatrix<T>& m2)
  { 
    TMVAssert(m1.rowsize() == m2.size());
    MultMM<false>(T(1),m1,BandMatrixViewOf(m2),m1); 
    return m1; 
  }

  template <class T> inline const BandMatrixView<CT>& operator*=(
      const BandMatrixView<CT>& m1, const GenDiagMatrix<T>& m2)
  {
    TMVAssert(m1.rowsize() == m2.size());
    MultMM<false>(T(1),m1,BandMatrixViewOf(m2),m1); 
    return m1; 
  }

  template <class T, class T1, class T2> 
    inline const BandMatrixView<T>& operator+=(
	const BandMatrixView<T>& m, const ProdDB<T,T1,T2>& pmm)
    { 
      TMVAssert(m.colsize() == pmm.colsize());
      TMVAssert(m.rowsize() == pmm.rowsize());
      TMVAssert(m.nhi() >= pmm.nhi());
      TMVAssert(m.nlo() >= pmm.nlo());
      MultMM<true>(pmm.GetX(),BandMatrixViewOf(pmm.GetM1()),pmm.GetM2(),m); 
      return m; 
    }

  template <class T> inline const BandMatrixView<CT>& operator+=(
	const BandMatrixView<CT>& m, const ProdDB<T,T,T>& pmm)
    { 
      TMVAssert(m.colsize() == pmm.colsize());
      TMVAssert(m.rowsize() == pmm.rowsize());
      TMVAssert(m.nhi() >= pmm.nhi());
      TMVAssert(m.nlo() >= pmm.nlo());
      MultMM<true>(pmm.GetX(),BandMatrixViewOf(pmm.GetM1()),pmm.GetM2(),m); 
      return m; 
    }

  template <class T, class T1, class T2> 
    inline const BandMatrixView<T>& operator-=(
	const BandMatrixView<T>& m, const ProdDB<T,T1,T2>& pmm)
    { 
      TMVAssert(m.colsize() == pmm.colsize());
      TMVAssert(m.rowsize() == pmm.rowsize());
      TMVAssert(m.nhi() >= pmm.nhi());
      TMVAssert(m.nlo() >= pmm.nlo());
      MultMM<true>(-pmm.GetX(),BandMatrixViewOf(pmm.GetM1()),pmm.GetM2(),m); 
      return m; 
    }

  template <class T> inline const BandMatrixView<CT>& operator-=(
	const BandMatrixView<CT>& m, const ProdDB<T,T,T>& pmm)
    { 
      TMVAssert(m.colsize() == pmm.colsize());
      TMVAssert(m.rowsize() == pmm.rowsize());
      TMVAssert(m.nhi() >= pmm.nhi());
      TMVAssert(m.nlo() >= pmm.nlo());
      MultMM<true>(-pmm.GetX(),BandMatrixViewOf(pmm.GetM1()),pmm.GetM2(),
	  m); 
      return m; 
    }

  template <class T, class T1, class T2> 
    inline const BandMatrixView<T>& operator+=(
	const BandMatrixView<T>& m, const ProdBD<T,T1,T2>& pmm)
    { 
      TMVAssert(m.colsize() == pmm.colsize());
      TMVAssert(m.rowsize() == pmm.rowsize());
      TMVAssert(m.nhi() >= pmm.nhi());
      TMVAssert(m.nlo() >= pmm.nlo());
      MultMM<true>(pmm.GetX(),pmm.GetM1(),BandMatrixViewOf(pmm.GetM2()),m); 
      return m; 
    }

  template <class T> inline const BandMatrixView<CT>& operator+=(
	const BandMatrixView<CT>& m, const ProdBD<T,T,T>& pmm)
    { 
      TMVAssert(m.colsize() == pmm.colsize());
      TMVAssert(m.rowsize() == pmm.rowsize());
      TMVAssert(m.nhi() >= pmm.nhi());
      TMVAssert(m.nlo() >= pmm.nlo());
      MultMM<true>(pmm.GetX(),pmm.GetM1(),BandMatrixViewOf(pmm.GetM2()),m); 
      return m; 
    }

  template <class T, class T1, class T2> 
    inline const BandMatrixView<T>& operator-=(
	const BandMatrixView<T>& m, const ProdBD<T,T1,T2>& pmm)
    {
      TMVAssert(m.colsize() == pmm.colsize());
      TMVAssert(m.rowsize() == pmm.rowsize());
      TMVAssert(m.nhi() >= pmm.nhi());
      TMVAssert(m.nlo() >= pmm.nlo());
      MultMM<true>(-pmm.GetX(),pmm.GetM1(),BandMatrixViewOf(pmm.GetM2()),m); 
      return m; 
    }

  template <class T> inline const BandMatrixView<CT>& operator-=(
	const BandMatrixView<CT>& m, const ProdBD<T,T,T>& pmm)
    {
      TMVAssert(m.colsize() == pmm.colsize());
      TMVAssert(m.rowsize() == pmm.rowsize());
      TMVAssert(m.nhi() >= pmm.nhi());
      TMVAssert(m.nlo() >= pmm.nlo());
      MultMM<true>(-pmm.GetX(),pmm.GetM1(),BandMatrixViewOf(pmm.GetM2()),
	  m);
      return m; 
    }

  template <class T> inline const DiagMatrixView<T>& operator*=(
      const DiagMatrixView<T>& m1, const GenBandMatrix<T>& m2)
  { 
    TMVAssert(m2.rowsize() == m1.size());
    TMVAssert(m2.colsize() == m1.size());
    TMVAssert(m2.nlo() == 0);
    TMVAssert(m2.nhi() == 0);
    MultMM<false>(T(1),BandMatrixViewOf(m1),m2,BandMatrixViewOf(m1)); 
    return m1; 
  }

  template <class T> inline const DiagMatrixView<CT>& operator*=(
      const DiagMatrixView<CT>& m1, const GenBandMatrix<T>& m2)
  {
    TMVAssert(m2.rowsize() == m1.size());
    TMVAssert(m2.colsize() == m1.size());
    TMVAssert(m2.nlo() == 0);
    TMVAssert(m2.nhi() == 0);
    MultMM<false>(T(1),BandMatrixViewOf(m1),m2,BandMatrixViewOf(m1)); 
    return m1; 
  }

  template <class T, class T1, class T2> 
    inline const DiagMatrixView<T>& operator+=(
	const DiagMatrixView<T>& m, const ProdDB<T,T1,T2>& pmm)
    { 
      TMVAssert(m.size() == pmm.colsize());
      TMVAssert(m.size() == pmm.rowsize());
      TMVAssert(pmm.nhi() == 0);
      TMVAssert(pmm.nlo() == 0);
      MultMM<true>(pmm.GetX(),BandMatrixViewOf(pmm.GetM1()),pmm.GetM2(),
	  BandMatrixViewOf(m)); 
      return m; 
    }

  template <class T> inline const DiagMatrixView<CT>& operator+=(
	const DiagMatrixView<CT>& m, const ProdDB<T,T,T>& pmm)
    { 
      TMVAssert(m.size() == pmm.colsize());
      TMVAssert(m.size() == pmm.rowsize());
      TMVAssert(pmm.nhi() == 0);
      TMVAssert(pmm.nlo() == 0);
      MultMM<true>(T(pmm.GetX()),BandMatrixViewOf(pmm.GetM1()),pmm.GetM2(),
	  BandMatrixViewOf(m)); 
      return m; 
    }

  template <class T, class T1, class T2> 
    inline const DiagMatrixView<T>& operator-=(
	const DiagMatrixView<T>& m, const ProdDB<T,T1,T2>& pmm)
    { 
      TMVAssert(m.size() == pmm.colsize());
      TMVAssert(m.size() == pmm.rowsize());
      TMVAssert(pmm.nhi() == 0);
      TMVAssert(pmm.nlo() == 0);
      MultMM<true>(-pmm.GetX(),BandMatrixViewOf(pmm.GetM1()),pmm.GetM2(),
	  BandMatrixViewOf(m)); 
      return m; 
    }

  template <class T> inline const DiagMatrixView<CT>& operator-=(
	const DiagMatrixView<CT>& m, const ProdDB<T,T,T>& pmm)
    { 
      TMVAssert(m.size() == pmm.colsize());
      TMVAssert(m.size() == pmm.rowsize());
      TMVAssert(pmm.nhi() == 0);
      TMVAssert(pmm.nlo() == 0);
      MultMM<true>(-pmm.GetX(),BandMatrixViewOf(pmm.GetM1()),pmm.GetM2(),
	  BandMatrixViewOf(m)); 
      return m; 
    }

  template <class T, class T1, class T2> 
    inline const DiagMatrixView<T>& operator+=(
	const DiagMatrixView<T>& m, const ProdBD<T,T1,T2>& pmm)
    { 
      TMVAssert(m.size() == pmm.colsize());
      TMVAssert(m.size() == pmm.rowsize());
      TMVAssert(pmm.nhi() == 0);
      TMVAssert(pmm.nlo() == 0);
      MultMM<true>(pmm.GetX(),pmm.GetM1(),BandMatrixViewOf(pmm.GetM2()),
	  BandMatrixViewOf(m)); 
      return m; 
    }

  template <class T> inline const DiagMatrixView<CT>& operator+=(
	const DiagMatrixView<CT>& m, const ProdBD<T,T,T>& pmm)
    { 
      TMVAssert(m.size() == pmm.colsize());
      TMVAssert(m.size() == pmm.rowsize());
      TMVAssert(pmm.nhi() == 0);
      TMVAssert(pmm.nlo() == 0);
      MultMM<true>(pmm.GetX(),pmm.GetM1(),BandMatrixViewOf(pmm.GetM2()),
	  BandMatrixViewOf(m)); 
      return m; 
    }

  template <class T, class T1, class T2> 
    inline const DiagMatrixView<T>& operator-=(
	const DiagMatrixView<T>& m, const ProdBD<T,T1,T2>& pmm)
    {
      TMVAssert(m.size() == pmm.colsize());
      TMVAssert(m.size() == pmm.rowsize());
      TMVAssert(pmm.nhi() == 0);
      TMVAssert(pmm.nlo() == 0);
      MultMM<true>(-pmm.GetX(),pmm.GetM1(),BandMatrixViewOf(pmm.GetM2()),
	  BandMatrixViewOf(m)); 
      return m; 
    }

  template <class T> inline const DiagMatrixView<CT>& operator-=(
	const DiagMatrixView<CT>& m, const ProdBD<T,T,T>& pmm)
    {
      TMVAssert(m.size() == pmm.colsize());
      TMVAssert(m.size() == pmm.rowsize());
      TMVAssert(pmm.nhi() == 0);
      TMVAssert(pmm.nlo() == 0);
      MultMM<true>(-pmm.GetX(),pmm.GetM1(),BandMatrixViewOf(pmm.GetM2()),
	  BandMatrixViewOf(m));
      return m; 
    }

#define PRODMM ProdBD
#define GENMATRIX1 GenBandMatrix
#define GENMATRIX2 GenDiagMatrix
#define PRODXM1 ProdXB
#define PRODXM2 ProdXD
#include "TMV_AuxProdMM.h"
#include "TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
  
#define PRODMM ProdDB
#define GENMATRIX1 GenDiagMatrix
#define GENMATRIX2 GenBandMatrix
#define PRODXM1 ProdXD
#define PRODXM2 ProdXB
#include "TMV_AuxProdMM.h"
#include "TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
  
  //
  // BandMatrix / % DiagMatrix
  //

  template <class T, class T1, class T2> class QuotBD : 
    public BandMatrixComposite<T>
  {
    public:
      inline QuotBD(const T _x, const GenBandMatrix<T1>& _m1,
	  const GenDiagMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert( m1.colsize() == m2.size() ); }
      inline size_t colsize() const { return m1.colsize(); }
      inline size_t rowsize() const { return m1.rowsize(); }
      inline int nlo() const { return m1.nlo(); }
      inline int nhi() const { return m1.nhi(); }
      inline StorageType stor() const { return BaseStorOf(m1); }
      inline T GetX() const { return x; }
      inline const GenBandMatrix<T1>& GetM1() const { return m1; }
      inline const GenDiagMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToB(const BandMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	TMVAssert(m0.nlo() >= nlo());
	TMVAssert(m0.nhi() >= nhi());
	MultMM<false>(x,BandMatrixViewOf(DiagMatrix<T2>(m2.Inverse())),m1,m0);
      }
      inline void AssignToB(const BandMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	TMVAssert(m0.nlo() >= nlo());
	TMVAssert(m0.nhi() >= nhi());
	MultMM<false>(x,BandMatrixViewOf(DiagMatrix<T2>(m2.Inverse())),m1,m0);
      }
    protected:
      const T x;
      const GenBandMatrix<T1>& m1;
      const GenDiagMatrix<T2>& m2;
  };

  template <class T, class T1, class T2> class RQuotBD : 
    public BandMatrixComposite<T>
  {
    public:
      inline RQuotBD(const T _x, const GenBandMatrix<T1>& _m1,
	  const GenDiagMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert( m1.rowsize() == m2.size() ); }
      inline size_t colsize() const { return m1.colsize(); }
      inline size_t rowsize() const { return m1.rowsize(); }
      inline int nlo() const { return m1.nlo(); }
      inline int nhi() const { return m1.nhi(); }
      inline StorageType stor() const { return BaseStorOf(m1); }
      inline T GetX() const { return x; }
      inline const GenBandMatrix<T1>& GetM1() const { return m1; }
      inline const GenDiagMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToB(const BandMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	TMVAssert(m0.nlo() >= nlo());
	TMVAssert(m0.nhi() >= nhi());
	MultMM<false>(x,m1,BandMatrixViewOf(DiagMatrix<T2>(m2.Inverse())),m0);
      }
      inline void AssignToB(const BandMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	TMVAssert(m0.nlo() >= nlo());
	TMVAssert(m0.nhi() >= nhi());
	MultMM<false>(x,m1,BandMatrixViewOf(DiagMatrix<T2>(m2.Inverse())),m0);
      }
    protected:
      const T x;
      const GenBandMatrix<T1>& m1;
      const GenDiagMatrix<T2>& m2;
  };

  template <class T> inline const BandMatrixView<T>& operator/=(
      const BandMatrixView<T>& m1, const GenDiagMatrix<T>& m2)
  { 
    TMVAssert(m1.colsize() == m2.size());
    MultMM<false>(T(1),BandMatrixViewOf(DiagMatrix<T>(m2.Inverse())),m1,m1); 
    return m1;
  }

  template <class T> inline const BandMatrixView<CT>& operator/=(
      const BandMatrixView<CT>& m1, const GenDiagMatrix<T>& m2)
  { 
    TMVAssert(m1.colsize() == m2.size());
    MultMM<false>(T(1),BandMatrixViewOf(DiagMatrix<T>(m2.Inverse())),m1,m1); 
    return m1;
  }

  template <class T> inline const BandMatrixView<T>& operator%=(
      const BandMatrixView<T>& m1, const GenDiagMatrix<T>& m2)
  {
    TMVAssert(m1.rowsize() == m2.size());
    MultMM<false>(T(1),m1,BandMatrixViewOf(DiagMatrix<T>(m2.Inverse())),m1); 
    return m1;
  }

  template <class T> inline const BandMatrixView<CT>& operator%=(
      const BandMatrixView<CT>& m1, const GenDiagMatrix<T>& m2)
  {
    TMVAssert(m1.rowsize() == m2.size());
    MultMM<false>(T(1),m1,BandMatrixViewOf(DiagMatrix<T>(m2.Inverse())),m1); 
    return m1;
  }

#define QUOTMM QuotBD
#define QUOTXM QuotXD
#define RQUOTMM RQuotBD
#define GENMATRIX1 GenBandMatrix
#define GENMATRIX2 GenDiagMatrix
#define PRODXM1 ProdXB
#define PRODXM2 ProdXD
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
#define GENMATRIX2 GenBandMatrix
#define PRODXM1 ProdXD
#define PRODXM2 ProdXB
#define QUOTXM QuotXB
#define TQUOTMM TransientQuotMB
#define TRQUOTMM TransientRQuotMB

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
