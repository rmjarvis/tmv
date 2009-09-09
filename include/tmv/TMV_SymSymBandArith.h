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


#ifndef TMV_SymSymBandArith_H
#define TMV_SymSymBandArith_H

#define CT std::complex<T>

namespace tmv {

  //
  // SymBandMatrix + SymMatrix
  //

  template <class T, class T1, class T2> 
  class SumsBS : 
    public SymMatrixComposite<T> 
  {
  public:
    inline SumsBS(T _x1, const GenSymBandMatrix<T1>& _m1, 
        T _x2, const GenSymMatrix<T2>& _m2) :
      x1(_x1),m1(_m1),x2(_x2),m2(_m2)
    { TMVAssert(m1.size() == m2.size()); }
    inline size_t size() const { return m2.size(); }
    inline SymType sym() const
    { return IsReal(T1()) ? m2.sym() : m1.sym(); }
    inline StorageType stor() const { return BaseStorOf(m2); }
    inline T GetX1() const { return x1; }
    inline const GenSymBandMatrix<T1>& GetM1() const { return m1; }
    inline T GetX2() const { return x2; }
    inline const GenSymMatrix<T2>& GetM2() const { return m2; }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    { 
      TMVAssert(IsReal(T()));
      TMVAssert(m0.colsize() == size());
      TMVAssert(m0.rowsize() == size());
      if (SameStorage(m0,m1)) {
        SymBandMatrix<T1> m1x = m1;
        MultXM(x2,m0=m2);
        AddMM(x1,m1x,m0);
      } else {
        MultXM(x2,m0=m2);
        AddMM(x1,m1,m0);
      }
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == size());
      TMVAssert(m0.rowsize() == size());
      if (SameStorage(m0,m1)) {
        if (m1.issym()) {
          SymBandMatrix<T1> m1x = m1;
          MultXM(x2,m0=m2);
          AddMM(x1,m1x,m0);
        } else {
          HermBandMatrix<T1> m1x = m1;
          MultXM(x2,m0=m2);
          AddMM(x1,m1x,m0);
        }
      } else {
        MultXM(x2,m0=m2);
        AddMM(x1,m1,m0);
      }
    }
    inline void AssignToS(const SymMatrixView<RealType(T)>& m0) const
    {
      TMVAssert(IsReal(T()));
      TMVAssert(m0.size() == size());
      TMVAssert(IsReal(T1()) || m0.issym() == m1.issym());
      TMVAssert(IsReal(T2()) || m0.issym() == m2.issym());
      TMVAssert(m0.issym() || IMAG(x1) == RealType(T)(0));
      TMVAssert(m0.issym() || IMAG(x2) == RealType(T)(0));
      if (SameStorage(m0,m1)) {
        SymBandMatrix<T1> m1x = m1;
        MultXM(x2,m0=m2);
        AddMM(x1,m1x,SymBandMatrixViewOf(m0,m1.nlo()));
      } else {
        MultXM(x2,m0=m2);
        AddMM(x1,m1,SymBandMatrixViewOf(m0,m1.nlo()));
      }
    }
    inline void AssignToS(const SymMatrixView<ComplexType(T)>& m0) const
    {
      TMVAssert(IsReal(T1()) || IsReal(T2()) || m1.sym() == m2.sym());
      TMVAssert(m0.size() == size());
      TMVAssert(IsReal(T1()) || m0.issym() == m1.issym());
      TMVAssert(IsReal(T2()) || m0.issym() == m2.issym());
      TMVAssert(m0.issym() || IMAG(x1) == RealType(T)(0));
      TMVAssert(m0.issym() || IMAG(x2) == RealType(T)(0));
      if (SameStorage(m0,m1)) {
        if (m1.issym()) {
          SymBandMatrix<T1> m1x = m1;
          MultXM(x2,m0=m2);
          AddMM(x1,m1x,SymBandMatrixViewOf(m0,m1.nlo()));
        } else {
          HermBandMatrix<T1> m1x = m1;
          MultXM(x2,m0=m2);
          AddMM(x1,m1x,SymBandMatrixViewOf(m0,m1.nlo()));
        }
      } else {
        MultXM(x2,m0=m2);
        AddMM(x1,m1,SymBandMatrixViewOf(m0,m1.nlo()));
      }
    }

  private:
    const T x1;
    const GenSymBandMatrix<T1>& m1;
    const T x2;
    const GenSymMatrix<T2>& m2;
  };

  template <class T> 
  inline const SymMatrixView<T>& operator+=(
      const SymMatrixView<T>& m1, const GenSymBandMatrix<T>& m2)
  {
    TMVAssert(m1.size() == m2.size());
    TMVAssert(IsReal(T()) || m1.sym() == m2.sym());
    AddMM(T(1),m2,SymBandMatrixViewOf(m1,m2.nlo()));
    return m1;
  }

  template <class T> 
  inline const SymMatrixView<CT>& operator+=(
      const SymMatrixView<CT>& m1, const GenSymBandMatrix<T>& m2)
  {
    TMVAssert(m1.size() == m2.size());
    AddMM(T(1),m2,SymBandMatrixViewOf(m1,m2.nlo()));
    return m1;
  }

  template <class T> 
  inline const SymMatrixView<T>& operator-=(
      const SymMatrixView<T>& m1, const GenSymBandMatrix<T>& m2)
  {
    TMVAssert(m1.size() == m2.size());
    TMVAssert(IsReal(T()) || m1.sym() == m2.sym());
    AddMM(T(-1),m2,SymBandMatrixViewOf(m1,m2.nlo()));
    return m1;
  }

  template <class T> 
  inline const SymMatrixView<CT>& operator-=(
      const SymMatrixView<CT>& m1, const GenSymBandMatrix<T>& m2)
  {
    TMVAssert(m1.size() == m2.size());
    AddMM(T(-1),m2,SymBandMatrixViewOf(m1,m2.nlo()));
    return m1;
  }

  template <class T, class T2> 
  inline const SymMatrixView<T>& operator+=(
      const SymMatrixView<T>& m, const ProdXsB<T,T2>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    TMVAssert(IsReal(T2()) || m.issym() == pxm.GetM().issym());
    TMVAssert(m.issym() || IMAG(pxm.GetX()) == RealType(T)(0));
    AddMM(pxm.GetX(),pxm.GetM(),SymBandMatrixViewOf(m,pxm.nlo()));
    return m;
  }

  template <class T> 
  inline const SymMatrixView<CT>& operator+=(
      const SymMatrixView<CT>& m, const ProdXsB<T,T>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    AddMM(pxm.GetX(),pxm.GetM(),SymBandMatrixViewOf(m,pxm.nlo()));
    return m;
  }

  template <class T, class T2> 
  inline const SymMatrixView<T>& operator-=(
      const SymMatrixView<T>& m, const ProdXsB<T,T2>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    TMVAssert(IsReal(T2()) || m.issym() == pxm.GetM().issym());
    TMVAssert(m.issym() || IMAG(pxm.GetX()) == RealType(T)(0));
    AddMM(-pxm.GetX(),pxm.GetM(),SymBandMatrixViewOf(m,pxm.nlo()));
    return m;
  }

  template <class T> 
  inline const SymMatrixView<CT>& operator-=(
      const SymMatrixView<CT>& m, const ProdXsB<T,T>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    AddMM(-pxm.GetX(),pxm.GetM(),SymBandMatrixViewOf(m,pxm.nlo()));
    return m;
  }


#define SUMMM SumsBS
#define GENMATRIX1 GenSymBandMatrix
#define GENMATRIX2 GenSymMatrix
#define PRODXM1 ProdXsB
#define PRODXM2 ProdXS
#include "tmv/TMV_AuxSumMM.h"
#include "tmv/TMV_AuxSumMMb.h"
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

  //
  // SymBandMatrix * SymMatrix
  //

  template <class T, class T1, class T2> 
  class ProdsBS : 
    public MatrixComposite<T>
  {
  public:
    inline ProdsBS(T _x, const GenSymBandMatrix<T1>& _m1,
        const GenSymMatrix<T2>& _m2) :
      x(_x), m1(_m1), m2(_m2)
    { TMVAssert(m1.size() == m2.size()); }
    inline size_t colsize() const { return m2.size(); }
    inline size_t rowsize() const { return m2.size(); }
    inline StorageType stor() const { return BaseStorOf(m2); }
    inline T GetX() const { return x; }
    inline const GenSymBandMatrix<T1>& GetM1() const { return m1; }
    inline const GenSymMatrix<T2>& GetM2() const { return m2; }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    {
      TMVAssert(IsReal(T()));
      TMVAssert(m0.colsize() == colsize());
      TMVAssert(m0.rowsize() == rowsize());
      Matrix<T> m2x = m2;
      MultMM<false>(x,m1,m2x,m0);
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    {
      TMVAssert(m0.colsize() == colsize());
      TMVAssert(m0.rowsize() == rowsize());
      Matrix<T> m2x = m2;
      MultMM<false>(x,m1,m2x,m0);
    }
  protected:
    T x;
    const GenSymBandMatrix<T1>& m1;
    const GenSymMatrix<T2>& m2;
  };

  template <class T, class T1, class T2> 
  class ProdSsB : 
    public MatrixComposite<T>
  {
  public:
    inline ProdSsB(T _x, const GenSymMatrix<T1>& _m1,
        const GenSymBandMatrix<T2>& _m2) :
      x(_x), m1(_m1), m2(_m2)
    { 
      TMVAssert(m2.colsize() == m1.size()); 
      TMVAssert(m2.rowsize() == m1.size()); 
    }
    inline size_t colsize() const { return m1.size(); }
    inline size_t rowsize() const { return m1.size(); }
    inline StorageType stor() const { return BaseStorOf(m1); }
    inline T GetX() const { return x; }
    inline const GenSymMatrix<T1>& GetM1() const { return m1; }
    inline const GenSymBandMatrix<T2>& GetM2() const { return m2; }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    {
      TMVAssert(IsReal(T()));
      TMVAssert(m0.colsize() == colsize());
      TMVAssert(m0.rowsize() == rowsize());
      Matrix<T> m1x = m1;
      MultMM<false>(x,m1x,m2,m0);
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    {
      TMVAssert(m0.colsize() == colsize());
      TMVAssert(m0.rowsize() == rowsize());
      Matrix<T> m1x = m1;
      MultMM<false>(x,m1x,m2,m0);
    }
  protected:
    T x;
    const GenSymMatrix<T1>& m1;
    const GenSymBandMatrix<T2>& m2;
  };

  template <class T, class T1, class T2>
  inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m, const ProdSsB<T,T1,T2>& pmm)
  {
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    AddMM(T(1),Matrix<T>(pmm),m);
    return m;
  }

  template <class T> 
  inline const MatrixView<CT>& operator+=(
      const MatrixView<CT>& m, const ProdSsB<T,T,T>& pmm)
  {
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    AddMM(T(1),Matrix<T>(pmm),m);
    return m;
  }

  template <class T, class T1, class T2>
  inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m, const ProdSsB<T,T1,T2>& pmm)
  {
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    AddMM(T(-1),Matrix<T>(pmm),m);
    return m;
  }

  template <class T> 
  inline const MatrixView<CT>& operator-=(
      const MatrixView<CT>& m, const ProdSsB<T,T,T>& pmm)
  {
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    AddMM(T(-1),Matrix<T>(pmm),m);
    return m;
  }

  template <class T, class T1, class T2>
  inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m, const ProdsBS<T,T1,T2>& pmm)
  {
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    AddMM(T(1),Matrix<T>(pmm),m);
    return m;
  }

  template <class T> 
  inline const MatrixView<CT>& operator+=(
      const MatrixView<CT>& m, const ProdsBS<T,T,T>& pmm)
  {
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    AddMM(T(1),Matrix<T>(pmm),m);
    return m;
  }

  template <class T, class T1, class T2>
  inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m, const ProdsBS<T,T1,T2>& pmm)
  {
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    AddMM(T(-1),Matrix<T>(pmm),m);
    return m;
  }

  template <class T> 
  inline const MatrixView<CT>& operator-=(
      const MatrixView<CT>& m, const ProdsBS<T,T,T>& pmm)
  {
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    AddMM(T(-1),Matrix<T>(pmm),m);
    return m;
  }


#define PRODMM ProdSsB
#define GENMATRIX1 GenSymMatrix
#define GENMATRIX2 GenSymBandMatrix
#define PRODXM1 ProdXS
#define PRODXM2 ProdXsB
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

#define PRODMM ProdsBS
#define GENMATRIX1 GenSymBandMatrix
#define GENMATRIX2 GenSymMatrix
#define PRODXM1 ProdXsB
#define PRODXM2 ProdXS
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

  //
  // SymMatrix / % SymBandMatrix
  //

#define GENMATRIX1 GenSymBandMatrix
#define GENMATRIX2 GenSymMatrix
#define PRODXM1 ProdXsB
#define PRODXM2 ProdXS
#define QUOTXM QuotXS
#define TQUOTMM TransientQuotMS
#define TRQUOTMM TransientRQuotMS
#include "tmv/TMV_AuxTQuotMM.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTXM
#undef TQUOTMM
#undef TRQUOTMM

#define GENMATRIX1 GenSymMatrix
#define GENMATRIX2 GenSymBandMatrix
#define PRODXM1 ProdXS
#define PRODXM2 ProdXsB
#define QUOTXM QuotXsB
#define TQUOTMM TransientQuotMsB
#define TRQUOTMM TransientRQuotMsB
#include "tmv/TMV_AuxTQuotMM.h"
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
