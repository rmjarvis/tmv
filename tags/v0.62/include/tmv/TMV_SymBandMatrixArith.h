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


#ifndef TMV_SymBandMatrixArith_H
#define TMV_SymBandMatrixArith_H

#include "tmv/TMV_BaseSymBandMatrix.h"
#include "tmv/TMV_SymBandMatrixArithFunc.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"

#define CT std::complex<T>
#define CCT ConjRef<std::complex<T> >
#define VCT VarConjRef<std::complex<T> >

namespace tmv {

  template <class T, class Tv> 
  class ProdXV;

  template <class T, class Tv> 
  class ProdXM;

  template <class T, UpLoType U, StorageType S, class Tx>
  inline SymBandMatrix<T,U,S>& operator+=(
      SymBandMatrix<T,U,S>& m, const Tx& x) 
  { m.View() += x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
  inline SymBandMatrix<T,U,S>& operator-=(
      SymBandMatrix<T,U,S>& m, const Tx& x) 
  { m.View() -= x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
  inline SymBandMatrix<T,U,S>& operator*=(
      SymBandMatrix<T,U,S>& m, const Tx& x) 
  { m.View() *= x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
  inline SymBandMatrix<T,U,S>& operator/=(
      SymBandMatrix<T,U,S>& m, const Tx& x) 
  { m.View() /= x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
  inline SymBandMatrix<T,U,S>& operator%=(
      SymBandMatrix<T,U,S>& m, const Tx& x) 
  { m.View() %= x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
  inline HermBandMatrix<T,U,S>& operator+=(
      HermBandMatrix<T,U,S>& m, const Tx& x) 
  { m.View() += x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
  inline HermBandMatrix<T,U,S>& operator-=(
      HermBandMatrix<T,U,S>& m, const Tx& x) 
  { m.View() -= x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
  inline HermBandMatrix<T,U,S>& operator*=(
      HermBandMatrix<T,U,S>& m, const Tx& x) 
  { m.View() *= x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
  inline HermBandMatrix<T,U,S>& operator/=(
      HermBandMatrix<T,U,S>& m, const Tx& x) 
  { m.View() /= x; return m; }

  template <class T, UpLoType U, StorageType S, class Tx>
  inline HermBandMatrix<T,U,S>& operator%=(
      HermBandMatrix<T,U,S>& m, const Tx& x) 
  { m.View() %= x; return m; }


  //
  // Scalar * SymBandMatrix
  //

  template <class T, class Tm> 
  class ProdXsB : 
    public SymBandMatrixComposite<T> 
  {
  public:
    inline ProdXsB(const T _x, const GenSymBandMatrix<Tm>& _m) :
      x(_x), m(_m) {}
    inline size_t size() const { return m.size(); }
    inline SymType sym() const { return m.issym() ? Sym : Herm; }
    inline int nlo() const { return m.nlo(); }
    inline StorageType stor() const { return BaseStorOf(m); }
    inline T GetX() const { return x; }
    inline const GenSymBandMatrix<Tm>& GetM() const { return m; }
    inline void AssignToB(const BandMatrixView<RealType(T)>& m0) const
    {
      TMVAssert(m0.colsize() == size());
      TMVAssert(m0.rowsize() == size());
      TMVAssert(m0.nlo() >= nlo());
      TMVAssert(m0.nhi() >= nlo());
      m0.Diags(0,m0.nhi()+1) = m.UpperBand();
      if (m.nlo() > 0) m0.Diags(-m0.nlo(),0) = m.LowerBandOff();
      else if (m0.nlo() > 0) m0.Diags(-m0.nlo(),0).Zero();
      MultXM(x,m0);
    }
    inline void AssignToB(const BandMatrixView<ComplexType(T)>& m0) const
    {
      TMVAssert(m0.colsize() == size());
      TMVAssert(m0.rowsize() == size());
      TMVAssert(m0.nlo() >= nlo());
      TMVAssert(m0.nhi() >= nlo());
      m0.Diags(0,m0.nhi()+1) = m.UpperBand();
      if (m.nlo() > 0) m0.Diags(-m0.nlo(),0) = m.LowerBandOff();
      else if (m0.nlo() > 0) m0.Diags(-m0.nlo(),0).Zero();
      MultXM(x,m0);
    }
    inline void AssignTosB(const SymBandMatrixView<RealType(T)>& m0) const
    {
      TMVAssert(m.issym() || IMAG(x)==RealType(T)(0) ); 
      TMVAssert(m0.size() == size());
      TMVAssert(m0.nlo() >= nlo());
      TMVAssert(IsReal(Tm()) || m0.issym() == m.issym());
      MultXM(x,m0=m);
    }
    inline void AssignTosB(const SymBandMatrixView<ComplexType(T)>& m0) const
    {
      TMVAssert(m.issym() || IMAG(x)==RealType(T)(0) ); 
      TMVAssert(m0.size() == size());
      TMVAssert(m0.nlo() >= nlo());
      TMVAssert(IsReal(Tm()) || m0.issym() == m.issym());
      MultXM(x,m0=m);
    }
  private:
    const T x;
    const GenSymBandMatrix<Tm>& m;
  };

  // m*=x
  template <class T> 
  inline const SymBandMatrixView<T>& operator*=(
      const SymBandMatrixView<T>& m, T x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == RealType(T)(0));
    MultXM(x,m);
    return m;
  }

  template <class T> 
  inline const SymBandMatrixView<CT>& operator*=(
      const SymBandMatrixView<CT>& m, T x) 
  { 
    MultXM(x,m);
    return m;
  }

  template <class T> 
  inline const SymBandMatrixView<CT>& operator*=(
      const SymBandMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == RealType(T)(0));
    MultXM(CT(x),m);
    return m;
  }

  template <class T> 
  inline const SymBandMatrixView<CT>& operator*=(
      const SymBandMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == RealType(T)(0));
    MultXM(CT(x),m);
    return m;
  }

  // m/=x
  template <class T> 
  inline const SymBandMatrixView<T>& operator/=(
      const SymBandMatrixView<T>& m, T x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == RealType(T)(0));
    MultXM(RealType(T)(1)/x,m);
    return m;
  }

  template <class T> 
  inline const SymBandMatrixView<CT>& operator/=(
      const SymBandMatrixView<CT>& m, T x) 
  { 
    MultXM(T(1)/x,m);
    return m;
  }

  template <class T> 
  inline const SymBandMatrixView<CT>& operator/=(
      const SymBandMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == RealType(T)(0));
    MultXM(T(1)/CT(x),m);
    return m;
  }

  template <class T> 
  inline const SymBandMatrixView<CT>& operator/=(
      const SymBandMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == RealType(T)(0));
    MultXM(T(1)/CT(x),m);
    return m;
  }

#define GENMATRIX GenSymBandMatrix
#define PRODXM ProdXsB
#include "tmv/TMV_AuxProdXM.h"
#undef GENMATRIX
#undef PRODXM

  //
  // SymBandMatrix + Scalar
  //

  template <class T, class Tm> 
  class SumsBX : 
    public SymBandMatrixComposite<T>
  {
  public:
    inline SumsBX(T _x1, const GenSymBandMatrix<Tm>& _m, T _x2) :
      x1(_x1), m(_m), x2(_x2) {}
    inline size_t size() const { return m.size(); }
    inline int nlo() const { return m.nlo(); }
    inline SymType sym() const { return m.issym() ? Sym : Herm; }
    inline StorageType stor() const { return BaseStorOf(m); }
    inline T GetX1() const { return x1; }
    inline const GenSymBandMatrix<Tm>& GetM() const { return m; }
    inline T GetX2() const { return x2; }
    inline void AssignToB(const BandMatrixView<RealType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == size());
      TMVAssert(m0.rowsize() == size());
      TMVAssert(m0.nlo() >= nlo());
      TMVAssert(m0.nhi() >= nlo());
      MultXM(x1,m0=m);
      m0.diag().AddToAll(REAL(x2));
    } 
    inline void AssignToB(const BandMatrixView<ComplexType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == size());
      TMVAssert(m0.rowsize() == size());
      TMVAssert(m0.nlo() >= nlo());
      TMVAssert(m0.nhi() >= nlo());
      MultXM(x1,m0=m);
      m0.diag().AddToAll(ComplexType(T)(x2));
    } 
    inline void AssignTosB(const SymBandMatrixView<RealType(T)>& m0) const
    { 
      TMVAssert(IsReal(Tm()) || m0.issym() == m.issym());
      TMVAssert(m.issym() || IMAG(x1)==RealType(T)(0));
      TMVAssert(m.issym() || IMAG(x2)==RealType(T)(0));
      TMVAssert(m0.size() == size());
      TMVAssert(m0.nlo() >= nlo());
      MultXM(x1,m0=m);
      m0.diag().AddToAll(REAL(x2));
    }
    inline void AssignTosB(const SymBandMatrixView<ComplexType(T)>& m0) const
    { 
      TMVAssert(IsReal(Tm()) || m0.issym() == m.issym());
      TMVAssert(m.issym() || IMAG(x1)==RealType(T)(0));
      TMVAssert(m.issym() || IMAG(x2)==RealType(T)(0));
      TMVAssert(m0.size() == size());
      TMVAssert(m0.nlo() >= nlo());
      MultXM(x1,m0=m);
      m0.diag().AddToAll(ComplexType(T)(x2));
    }
  private:
    const T x1;
    const GenSymBandMatrix<Tm>& m;
    const T x2;
  };

  // m+=x
  template <class T> 
  inline const SymBandMatrixView<T>& operator+=(
      const SymBandMatrixView<T>& m, T x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == RealType(T)(0));
    m.diag().AddToAll(x);
    return m; 
  }

  template <class T> 
  inline const SymBandMatrixView<CT>& operator+=(
      const SymBandMatrixView<CT>& m, T x) 
  { 
    m.diag().AddToAll(CT(x));
    return m; 
  }

  template <class T> 
  inline const SymBandMatrixView<CT>& operator+=(
      const SymBandMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == T(0));
    m.diag().AddToAll(CT(x));
    return m; 
  }

  template <class T> 
  inline const SymBandMatrixView<CT>& operator+=(
      const SymBandMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == T(0));
    m.diag().AddToAll(CT(x));
    return m; 
  }

  // m-=x
  template <class T> 
  inline const SymBandMatrixView<T>& operator-=(
      const SymBandMatrixView<T>& m, T x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == RealType(T)(0));
    m.diag().AddToAll(-x);
    return m; 
  }

  template <class T> 
  inline const SymBandMatrixView<CT>& operator-=(
      const SymBandMatrixView<CT>& m, T x) 
  { 
    m.diag().AddToAll(CT(-x));
    return m; 
  }

  template <class T> 
  inline const SymBandMatrixView<CT>& operator-=(
      const SymBandMatrixView<CT>& m, CCT x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == T(0));
    m.diag().AddToAll(-CT(x));
    return m; 
  }

  template <class T> 
  inline const SymBandMatrixView<CT>& operator-=(
      const SymBandMatrixView<CT>& m, VCT x) 
  { 
    TMVAssert(m.issym() || IMAG(x) == T(0));
    m.diag().AddToAll(-CT(x));
    return m; 
  }

#define SUMMX SumsBX
#define GENMATRIX GenSymBandMatrix
#define PRODXM ProdXsB
#include "tmv/TMV_AuxSumMX.h"
#undef SUMMX 
#undef GENMATRIX
#undef PRODXM

  //
  // SymBandMatrix + SymBandMatrix
  //

  template <class T, class T1, class T2> 
  class SumsBsB : 
    public SymBandMatrixComposite<T> 
  {
  public:
    inline SumsBsB(T _x1, const GenSymBandMatrix<T1>& _m1, 
        T _x2, const GenSymBandMatrix<T2>& _m2) :
      x1(_x1),m1(_m1),x2(_x2),m2(_m2)
    { TMVAssert(m1.size() == m2.size()); }
    inline size_t size() const { return m1.size(); }
    inline int nlo() const { return MAX(m1.nlo(),m2.nlo()); }
    inline SymType sym() const 
    { return IsReal(T1()) ? m2.sym() : m1.sym(); }
    inline StorageType stor() const 
    { return m1.stor() == m2.stor() ? BaseStorOf(m1) : DiagMajor; }
    inline T GetX1() const { return x1; }
    inline const GenSymBandMatrix<T1>& GetM1() const { return m1; }
    inline T GetX2() const { return x2; }
    inline const GenSymBandMatrix<T2>& GetM2() const { return m2; }
    inline void AssignToB(const BandMatrixView<RealType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == size());
      TMVAssert(m0.rowsize() == size());
      TMVAssert(m0.nlo() >= nlo());
      TMVAssert(m0.nhi() >= nlo());
      AddMM(x1,m1,x2,m2,m0);
    } 
    inline void AssignToB(const BandMatrixView<ComplexType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == size());
      TMVAssert(m0.rowsize() == size());
      TMVAssert(m0.nlo() >= nlo());
      TMVAssert(m0.nhi() >= nlo());
      AddMM(x1,m1,x2,m2,m0);
    } 
    inline void AssignTosB(const SymBandMatrixView<RealType(T)>& m0) const
    { 
      TMVAssert(IsReal(T1()) || IsReal(T2()) || m1.sym() == m2.sym());
      TMVAssert(m0.size() == size());
      TMVAssert(m0.nlo() >= nlo());
      TMVAssert(IsReal(T1()) || m0.issym() == m1.issym());
      TMVAssert(IsReal(T2()) || m0.issym() == m2.issym());
      TMVAssert(m0.issym() || IMAG(x1) == RealType(T)(0));
      TMVAssert(m0.issym() || IMAG(x2) == RealType(T)(0));
      AddMM(x1,m1,x2,m2,m0);
    }
    inline void AssignTosB(const SymBandMatrixView<ComplexType(T)>& m0) const
    { 
      TMVAssert(IsReal(T1()) || IsReal(T2()) || m1.sym() == m2.sym());
      TMVAssert(m0.size() == size());
      TMVAssert(m0.nlo() >= nlo());
      TMVAssert(IsReal(T1()) || m0.issym() == m1.issym());
      TMVAssert(IsReal(T2()) || m0.issym() == m2.issym());
      TMVAssert(m0.issym() || IMAG(x1) == RealType(T)(0));
      TMVAssert(m0.issym() || IMAG(x2) == RealType(T)(0));
      AddMM(x1,m1,x2,m2,m0);
    }
  private:
    T x1;
    const GenSymBandMatrix<T1>& m1;
    T x2;
    const GenSymBandMatrix<T2>& m2;
  };

  template <class T> 
  inline const SymBandMatrixView<T>& operator+=(
      const SymBandMatrixView<T>& m1, const GenSymBandMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(m1.nlo() >= m2.nlo());
    TMVAssert(IsReal(T()) || m1.sym() == m2.sym());
    AddMM(T(1),m2,m1); return m1; 
  }

  template <class T> 
  inline const SymBandMatrixView<CT>& operator+=(
      const SymBandMatrixView<CT>& m1, const GenSymBandMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(m1.nlo() >= m2.nlo());
    AddMM(T(1),m2,m1); return m1; 
  }

  template <class T> 
  inline const SymBandMatrixView<T>& operator-=(
      const SymBandMatrixView<T>& m1, const GenSymBandMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(m1.nlo() >= m2.nlo());
    TMVAssert(IsReal(T()) || m1.sym() == m2.sym());
    AddMM(T(-1),m2,m1); return m1; 
  }

  template <class T> 
  inline const SymBandMatrixView<CT>& operator-=(
      const SymBandMatrixView<CT>& m1, const GenSymBandMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(m1.nlo() >= m2.nlo());
    AddMM(T(-1),m2,m1); return m1; 
  }

  template <class T, class T2> 
  inline const SymBandMatrixView<T>& operator+=(
      const SymBandMatrixView<T>& m, const ProdXsB<T,T2>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    TMVAssert(m.nlo() >= pxm.nlo());
    TMVAssert(IsReal(T2()) || m.issym() == pxm.GetM().issym());
    TMVAssert(m.issym() || IMAG(pxm.GetX()) == RealType(T)(0));
    AddMM(pxm.GetX(),pxm.GetM(),m);
    return m;
  }

  template <class T> 
  inline const SymBandMatrixView<CT>& operator+=(
      const SymBandMatrixView<CT>& m, const ProdXsB<T,T>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    TMVAssert(m.nlo() >= pxm.nlo());
    AddMM(pxm.GetX(),pxm.GetM(),m);
    return m;
  }

  template <class T, class T2> 
  inline const SymBandMatrixView<T>& operator-=(
      const SymBandMatrixView<T>& m, const ProdXsB<T,T2>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    TMVAssert(m.nlo() >= pxm.nlo());
    TMVAssert(IsReal(T2()) || m.issym() == pxm.GetM().issym());
    TMVAssert(m.issym() || IMAG(pxm.GetX()) == RealType(T)(0));
    AddMM(-pxm.GetX(),pxm.GetM(),m);
    return m;
  }

  template <class T> 
  inline const SymBandMatrixView<CT>& operator-=(
      const SymBandMatrixView<CT>& m, const ProdXsB<T,T>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    TMVAssert(m.nlo() >= pxm.nlo());
    AddMM(-pxm.GetX(),pxm.GetM(),m);
    return m;
  }

#define SUMMM SumsBsB
#define GENMATRIX1 GenSymBandMatrix
#define GENMATRIX2 GenSymBandMatrix
#define PRODXM1 ProdXsB
#define PRODXM2 ProdXsB
#include "tmv/TMV_AuxSumMM.h"
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2


  //
  // SymBandMatrix * SymBandMatrix
  //

  template <class T, class T1, class T2> 
  class ProdsBsB : 
    public BandMatrixComposite<T>
  {
  public:
    inline ProdsBsB(T _x, const GenSymBandMatrix<T1>& _m1,
        const GenSymBandMatrix<T2>& _m2) :
      x(_x), m1(_m1), m2(_m2)
    { TMVAssert(m1.size() == m2.size()); }
    inline size_t colsize() const { return m1.size(); }
    inline size_t rowsize() const { return m1.size(); }
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
    inline const GenSymBandMatrix<T1>& GetM1() const { return m1; }
    inline const GenSymBandMatrix<T2>& GetM2() const { return m2; }
    inline void AssignToB(const BandMatrixView<RealType(T)>& m0) const
    {
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
    const GenSymBandMatrix<T1>& m1;
    const GenSymBandMatrix<T2>& m2;
  };

  template <class T, class T1, class T2> 
  inline const BandMatrixView<T>& operator+=(
      const BandMatrixView<T>& m, const ProdsBsB<T,T1,T2>& pmm)
  {
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    TMVAssert(m.nlo() >= pmm.nlo());;
    TMVAssert(m.nhi() >= pmm.nhi());;
    MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
    return m; 
  }

  template <class T> 
  inline const BandMatrixView<CT>& operator+=(
      const BandMatrixView<CT>& m, const ProdsBsB<T,T,T>& pmm)
  {
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    TMVAssert(m.nlo() >= pmm.nlo());;
    TMVAssert(m.nhi() >= pmm.nhi());;
    MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
    return m; 
  }

  template <class T, class T1, class T2> 
  inline const BandMatrixView<T>& operator-=(
      const BandMatrixView<T>& m, const ProdsBsB<T,T1,T2>& pmm)
  { 
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    TMVAssert(m.nlo() >= pmm.nlo());;
    TMVAssert(m.nhi() >= pmm.nhi());;
    MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
    return m; 
  }

  template <class T> 
  inline const BandMatrixView<CT>& operator-=(
      const BandMatrixView<CT>& m, const ProdsBsB<T,T,T>& pmm)
  { 
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    TMVAssert(m.nlo() >= pmm.nlo());;
    TMVAssert(m.nhi() >= pmm.nhi());;
    MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
    return m; 
  }

  template <class T, class T1, class T2> 
  inline const MatrixView<T>& operator+=(
      const MatrixView<T>& m, const ProdsBsB<T,T1,T2>& pmm)
  { 
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    BandMatrixViewOf(m,pmm.nlo(),pmm.nhi()) += pmm; 
    return m;
  }

  template <class T> 
  inline const MatrixView<CT>& operator+=(
      const MatrixView<CT>& m, const ProdsBsB<T,T,T>& pmm)
  { 
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    BandMatrixViewOf(m,pmm.nlo(),pmm.nhi()) += pmm; 
    return m;
  }

  template <class T, class T1, class T2> 
  inline const MatrixView<T>& operator-=(
      const MatrixView<T>& m, const ProdsBsB<T,T1,T2>& pmm)
  {
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    BandMatrixViewOf(m,pmm.nlo(),pmm.nhi()) -= pmm;
    return m;
  }

  template <class T> 
  inline const MatrixView<CT>& operator-=(
      const MatrixView<CT>& m, const ProdsBsB<T,T,T>& pmm)
  {
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    BandMatrixViewOf(m,pmm.nlo(),pmm.nhi()) -= pmm; 
    return m;
  }

#define PRODMM ProdsBsB
#define GENMATRIX1 GenSymBandMatrix
#define GENMATRIX2 GenSymBandMatrix
#define PRODXM1 ProdXsB
#define PRODXM2 ProdXsB
#include "tmv/TMV_AuxProdMM.h"
#include "tmv/TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2


#define GENMATRIX GenSymBandMatrix
#define PRODXM ProdXsB
#define PRODMV ProdsBV
#define PRODVM ProdVsB
#define QUOTVM QuotVsB
#define RQUOTVM RQuotVsB
#define QUOTXM QuotXsB
#include "tmv/TMV_AuxVecComposite.h"
#undef GENMATRIX
#undef PRODXM
#undef PRODMV
#undef PRODVM
#undef QUOTVM
#undef RQUOTVM
#undef QUOTXM


#define GENMATRIX GenSymBandMatrix
#define GENMATRIX1 GenMatrix
#define GENMATRIX2 GenSymBandMatrix
#define PRODXM1 ProdXM
#define PRODXM2 ProdXsB
#define SUMMM SumMsB
#define PRODMM ProdMsB
#define QUOTXM QuotXsB
#define QUOTMM QuotMsB
#define RQUOTMM RQuotMsB
#define TQUOTMM TransientQuotMsB
#define TRQUOTMM TransientRQuotMsB

#include "tmv/TMV_AuxMatComposite1.h"
#include "tmv/TMV_AuxSumMMb.h"
#include "tmv/TMV_AuxMatComposite3.h"

  // sB/sB -> Matrix.  Use TransientQuot
#undef GENMATRIX
#undef GENMATRIX1
#undef PRODXM1
#define GENMATRIX1 GenSymBandMatrix
#define PRODXM1 ProdXsB
#include "tmv/TMV_AuxTQuotMM.h"

#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef SUMMM
#undef PRODMM
#undef QUOTMM
#undef RQUOTMM
#undef TQUOTMM
#undef TRQUOTMM
#undef QUOTXM


#define GENMATRIX1 GenSymBandMatrix
#define GENMATRIX2 GenMatrix
#define PRODXM1 ProdXsB
#define PRODXM2 ProdXM
#define SUMMMa SumMsB
#define SUMMM SumsBM
#define PRODMM ProdsBM
#define TQUOTMM TransientQuotMM
#define TRQUOTMM TransientRQuotMM
#define QUOTXM QuotXM

#include "tmv/TMV_AuxMatComposite2.h"
#include "tmv/TMV_AuxTQuotMM.h"

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
#include "tmv/TMV_DiagSymBandArith.h"
#endif

#ifdef TMV_TriMatrixArith_H
#include "tmv/TMV_TriSymBandArith.h"
#endif

#ifdef TMV_BandMatrixArith_H
#include "tmv/TMV_BandSymBandArith.h"
#endif

#ifdef TMV_SymMatrixArith_H
#include "tmv/TMV_SymSymBandArith.h"
#endif

#endif
