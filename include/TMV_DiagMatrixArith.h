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


#ifndef TMV_DiagMatrixArith_H
#define TMV_DiagMatrixArith_H

#include "TMV_BaseDiagMatrix.h"
#include "TMV_DiagMatrixArithFunc.h"
#include "TMV_VectorArithFunc.h"
#include "TMV_MatrixArithFunc.h"

#define CT std::complex<T>
#define CCT ConjRef<std::complex<T> >
#define VCT VarConjRef<std::complex<T> >

namespace tmv {

  template <class T, class Tv> class ProdXV;
  template <class T, class Tv> class ProdXM;

  template <class T, IndexStyle I, class Tx> inline DiagMatrix<T,I>& operator+=(
      DiagMatrix<T,I>& m, const Tx& x)
  { m.View() += x; return m; }

  template <class T, IndexStyle I, class Tx> inline DiagMatrix<T,I>& operator-=(
      DiagMatrix<T,I>& m, const Tx& x)
  { m.View() -= x; return m; }

  template <class T, IndexStyle I, class Tx> inline DiagMatrix<T,I>& operator*=(
      DiagMatrix<T,I>& m, const Tx& x)
  { m.View() *= x; return m; }

  template <class T, IndexStyle I, class Tx> inline DiagMatrix<T,I>& operator/=(
      DiagMatrix<T,I>& m, const Tx& x)
  { m.View() /= x; return m; }

  template <class T, IndexStyle I, class Tx> inline DiagMatrix<T,I>& operator%=(
      DiagMatrix<T,I>& m, const Tx& x)
  { m.View() %= x; return m; }


  //
  // Scalar * DiagMatrix
  //

  template <class T, class T2> class ProdXD : 
    public DiagMatrixComposite<T> 
  {
    public:
      inline ProdXD(T _x, const GenDiagMatrix<T2>& _m) : x(_x), m(_m) {}
      inline size_t size() const { return m.size(); }
      inline T GetX() const { return x; }
      inline const GenDiagMatrix<T2>& GetM() const { return m; }
      inline void AssignToD(const DiagMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	MultXV(x,m.diag(),m0.diag());
      }
      inline void AssignToD(const DiagMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.size() == size());
	MultXV(x,m.diag(),m0.diag());
      }
    private:
      const T x;
      const GenDiagMatrix<T2>& m;
  };

  template <class T> inline const DiagMatrixView<T>& operator*=(
      const DiagMatrixView<T>& m, T x) 
  { MultXV(x,m.diag()); return m; }

  template <class T> inline const DiagMatrixView<CT>& operator*=(
      const DiagMatrixView<CT>& m, T x) 
  { MultXV(x,m.diag()); return m; }

  template <class T> inline const DiagMatrixView<CT>& operator*=(
      const DiagMatrixView<CT>& m, CCT x) 
  { MultXV(CT(x),m.diag()); return m; }

  template <class T> inline const DiagMatrixView<CT>& operator*=(
      const DiagMatrixView<CT>& m, VCT x) 
  { MultXV(CT(x),m.diag()); return m; }

  template <class T> inline const DiagMatrixView<T>& operator/=(
      const DiagMatrixView<T>& m, T x) 
  { MultXV(RealType(T)(1)/x,m.diag()); return m; }

  template <class T> inline const DiagMatrixView<CT>& operator/=(
      const DiagMatrixView<CT>& m, T x) 
  { MultXV(T(1)/x,m.diag()); return m; }

  template <class T> inline const DiagMatrixView<CT>& operator/=(
      const DiagMatrixView<CT>& m, CCT x) 
  { MultXV(T(1)/CT(x),m.diag()); return m; }

  template <class T> inline const DiagMatrixView<CT>& operator/=(
      const DiagMatrixView<CT>& m, VCT x) 
  { MultXV(T(1)/CT(x),m.diag()); return m; }

#define GENMATRIX GenDiagMatrix
#define PRODXM ProdXD
#include "TMV_AuxProdXM.h"
  // Defines things like -m, x*m, m*x, x*(x*m), etc.
#undef GENMATRIX
#undef PRODXM

  //
  // DiagMatrix + Scalar
  //

  template <class T, class T2> class SumDX : 
    public DiagMatrixComposite<T>
  {
    public:
      inline SumDX(T _x1, const GenDiagMatrix<T2>& _m, T _x2) :
	x1(_x1), m(_m), x2(_x2) {}
      inline size_t size() const { return m.size(); }
      inline T GetX1() const { return x1; }
      inline const GenDiagMatrix<T2>& GetM() const { return m; }
      inline T GetX2() const { return x2; }
      inline void AssignToD(const DiagMatrixView<RealType(T)>& m0) const
      { 
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	MultXV(x1,m.diag(),m0.diag());
	m0.diag().AddToAll(REAL(x2));
      }
      inline void AssignToD(const DiagMatrixView<ComplexType(T)>& m0) const
      { 
	TMVAssert(m0.size() == size());
	MultXV(x1,m.diag(),m0.diag());
	m0.diag().AddToAll(ComplexType(T)(x2));
      }
    private:
      const T x1;
      const GenDiagMatrix<T2>& m;
      const T x2;
  };

  template <class T> inline const DiagMatrixView<T>& operator+=(
      const DiagMatrixView<T>& m, T x) 
  { m.diag().AddToAll(x); return m; }

  template <class T> inline const DiagMatrixView<CT>& operator+=(
      const DiagMatrixView<CT>& m, T x) 
  { m.diag().AddToAll(CT(x)); return m; }

  template <class T> inline const DiagMatrixView<CT>& operator+=(
      const DiagMatrixView<CT>& m, CCT x) 
  { m.diag().AddToAll(CT(x)); return m; }

  template <class T> inline const DiagMatrixView<CT>& operator+=(
      const DiagMatrixView<CT>& m, VCT x) 
  { m.diag().AddToAll(CT(x)); return m; }

  template <class T> inline const DiagMatrixView<T>& operator-=(
      const DiagMatrixView<T>& m, T x) 
  { m.diag().AddToAll(-x); return m; }

  template <class T> inline const DiagMatrixView<CT>& operator-=(
      const DiagMatrixView<CT>& m, T x) 
  { m.diag().AddToAll(CT(-x)); return m; }

  template <class T> inline const DiagMatrixView<CT>& operator-=(
      const DiagMatrixView<CT>& m, CCT x) 
  { m.diag().AddToAll(-CT(x)); return m; }

  template <class T> inline const DiagMatrixView<CT>& operator-=(
      const DiagMatrixView<CT>& m, VCT x) 
  { m.diag().AddToAll(-CT(x)); return m; }

#define SUMMX SumDX
#define GENMATRIX GenDiagMatrix
#define PRODXM ProdXD
#include "TMV_AuxSumMX.h"
  // Defines things like m+x, x+m, x-m, m-x, x+x*m, x*(x+m), etc.
#undef SUMMX
#undef GENMATRIX
#undef PRODXM

  //
  // DiagMatrix + DiagMatrix
  //
  
  template <class T, class T1, class T2> class SumDD : 
    public DiagMatrixComposite<T> 
  {
    public:
      inline SumDD(T _x1, const GenDiagMatrix<T1>& _m1, 
	  T _x2, const GenDiagMatrix<T2>& _m2) :
	x1(_x1),m1(_m1),x2(_x2),m2(_m2)
      { TMVAssert(m1.size() == m2.size()); }
      inline size_t size() const { return m1.size(); }
      inline T GetX1() const { return x1; }
      inline const GenDiagMatrix<T1>& GetM1() const { return m1; }
      inline T GetX2() const { return x2; }
      inline const GenDiagMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToD(const DiagMatrixView<RealType(T)>& m0) const
      { 
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	AddVV(x1,m1.diag(),x2,m2.diag(),m0.diag());
      }
      inline void AssignToD(const DiagMatrixView<ComplexType(T)>& m0) const
      { 
	TMVAssert(m0.size() == size());
	AddVV(x1,m1.diag(),x2,m2.diag(),m0.diag());
      }
    private:
      T x1;
      const GenDiagMatrix<T1>& m1;
      T x2;
      const GenDiagMatrix<T2>& m2;
  };

  template <class T> inline const DiagMatrixView<T>& operator+=(
      const DiagMatrixView<T>& m1, const GenDiagMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    AddVV(T(1),m2.diag(),m1.diag());
    return m1; 
  }

  template <class T> inline const DiagMatrixView<CT>& operator+=(
      const DiagMatrixView<CT>& m1, const GenDiagMatrix<T>& m2) 
  {
    TMVAssert(m1.size() == m2.size());
    AddVV(T(1),m2.diag(),m1.diag());
    return m1; 
  }

  template <class T> inline const DiagMatrixView<T>& operator-=(
      const DiagMatrixView<T>& m1, const GenDiagMatrix<T>& m2) 
  { 
    TMVAssert(m1.size() == m2.size());
    AddVV(T(-1),m2.diag(),m1.diag());
    return m1; 
  }

  template <class T> inline const DiagMatrixView<CT>& operator-=(
      const DiagMatrixView<CT>& m1, const GenDiagMatrix<T>& m2) 
  {
    TMVAssert(m1.size() == m2.size());
    AddVV(T(-1),m2.diag(),m1.diag());
    return m1; 
  }

  template <class T, class T2> inline const DiagMatrixView<T>& operator+=(
      const DiagMatrixView<T>& m, const ProdXD<T,T2>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    AddVV(pxm.GetX(),pxm.GetM().diag(),m.diag());
    return m;
  }

  template <class T> inline const DiagMatrixView<CT>& operator+=(
      const DiagMatrixView<CT>& m, const ProdXD<T,T>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    AddVV(pxm.GetX(),pxm.GetM().diag(),m.diag());
    return m;
  }

  template <class T, class T2> inline const DiagMatrixView<T>& operator-=(
      const DiagMatrixView<T>& m, const ProdXD<T,T2>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    AddVV(-pxm.GetX(),pxm.GetM().diag(),m.diag());
    return m;
  }

  template <class T> inline const DiagMatrixView<CT>& operator-=(
      const DiagMatrixView<CT>& m, const ProdXD<T,T>& pxm)
  {
    TMVAssert(m.size() == pxm.size());
    AddVV(-pxm.GetX(),pxm.GetM().diag(),m.diag());
    return m;
  }

#define SUMMM SumDD
#define GENMATRIX1 GenDiagMatrix
#define GENMATRIX2 GenDiagMatrix
#define PRODXM1 ProdXD
#define PRODXM2 ProdXD
#include "TMV_AuxSumMM.h"
  // Defines things like m+m, m-m, -(m+m), x*(m+m), x*m+x*m, etc.
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

  //
  // DiagMatrix * DiagMatrix
  //

  template <class T, class T1, class T2> class ProdDD : 
    public DiagMatrixComposite<T>
  {
    public:
      inline ProdDD(T _x, const GenDiagMatrix<T1>& _m1,
	  const GenDiagMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert( m1.size() == m2.size()); }
      inline size_t size() const { return m1.size(); }
      inline T GetX() const { return x; }
      inline const GenDiagMatrix<T1>& GetM1() const { return m1; }
      inline const GenDiagMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToD(const DiagMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	MultMM<false>(x,m1,m2,m0);
      }
      inline void AssignToD(const DiagMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.size() == size());
	MultMM<false>(x,m1,m2,m0);
      }
    private:
      T x;
      const GenDiagMatrix<T1>& m1;
      const GenDiagMatrix<T2>& m2;
  };

  template <class T> inline const DiagMatrixView<T>& operator*=(
      const DiagMatrixView<T>& m1, const GenDiagMatrix<T>& m2)
  {
    TMVAssert(m1.size() == m2.size());
    MultMM<false>(T(1),m1,m2,m1); 
    return m1; 
  }

  template <class T> inline const DiagMatrixView<CT>& operator*=(
      const DiagMatrixView<CT>& m1, const GenDiagMatrix<T>& m2)
  {
    TMVAssert(m1.size() == m2.size());
    MultMM<false>(T(1),m1,m2,m1); 
    return m1; 
  }

  template <class T, class T2> inline const DiagMatrixView<T>& operator*=(
      const DiagMatrixView<T>& m, const ProdXD<T,T2>& pxm)
  { 
    TMVAssert(m.size()==pxm.size());
    MultMM<false>(pxm.GetX(),m,pxm.GetM(),m); 
    return m; 
  }

  template <class T> inline const DiagMatrixView<CT>& operator*=(
      const DiagMatrixView<CT>& m, const ProdXD<T,T>& pxm)
  { 
    TMVAssert(m.size()==pxm.size());
    MultMM<false>(pxm.GetX(),m,pxm.GetM(),m); 
    return m; 
  }

  template <class T, class T2, class T3> 
    inline const DiagMatrixView<T>& operator+=(
	const DiagMatrixView<T>& m, const ProdDD<T,T2,T3>& pmm)
    { 
      TMVAssert(m.size() == pmm.size());
      MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T> inline const DiagMatrixView<CT>& operator+=(
      const DiagMatrixView<CT>& m, const ProdDD<T,T,T>& pmm)
  { 
    TMVAssert(m.size() == pmm.size());
    MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
    return m; 
  }

  template <class T, class T2, class T3> 
    inline const DiagMatrixView<T>& operator-=(
	const DiagMatrixView<T>& m, const ProdDD<T,T2,T3>& pmm)
    { 
      TMVAssert(m.size() == pmm.size());
      MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
      return m; 
    }

  template <class T> inline const DiagMatrixView<CT>& operator-=(
      const DiagMatrixView<CT>& m, const ProdDD<T,T,T>& pmm)
  { 
    TMVAssert(m.colsize() == pmm.colsize());
    TMVAssert(m.rowsize() == pmm.rowsize());
    MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m); 
    return m; 
  }

#define PRODMM ProdDD
#define GENMATRIX1 GenDiagMatrix
#define GENMATRIX2 GenDiagMatrix
#define PRODXM1 ProdXD
#define PRODXM2 ProdXD
#include "TMV_AuxProdMM.h"
  // Defines things like m*m, m*(x*m), (x*m)*(x*m), etc.
#include "TMV_AuxProdMMa.h"
  // Defines things like (x*m*m)*x, -(x*m*m), etc.
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2


  //
  // Scalar / DiagMatrix
  //

  template <class T, class Tm> class QuotXD : 
    public DiagMatrixComposite<T>
  {
    public:
      inline QuotXD(T _x, const GenDiagMatrix<Tm>& _m) :
	x(_x), m(_m) {}
      inline size_t size() const { return m.size(); }
      inline T GetX() const { return x; }
      inline const GenDiagMatrix<Tm>& GetM() const { return m; }
      inline void AssignToD(const DiagMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	(m0=m).Inverse(m0);
	if (x != T(1)) MultXV(x,m0.diag());
      }
      inline void AssignToD(const DiagMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.size() == size());
	(m0=m).Inverse(m0);
	if (x != T(1)) MultXV(x,m0.diag());
      }
    protected:
      T x;
      const GenDiagMatrix<Tm>& m;
  };

#define QUOTXM QuotXD
#define GENMATRIX GenDiagMatrix
#define PRODXM ProdXD
#include "TMV_AuxQuotXM.h"
  // Defines x/m, x%m
#undef GENMATRIX
#undef PRODXM
  

  //
  // DiagMatrix / % DiagMatrix
  //

  template <class T, class T1, class T2> class QuotDD : 
    public DiagMatrixComposite<T>
  {
    public:
      inline QuotDD(const T _x, const GenDiagMatrix<T1>& _m1,
	  const GenDiagMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert( m1.size() == m2.size() ); }
      inline size_t size() const { return m1.size(); }
      inline T GetX() const { return x; }
      inline const GenDiagMatrix<T1>& GetM1() const { return m1; }
      inline const GenDiagMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToD(const DiagMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.size() == size());
	m2.Div(m1,m0);
	if (x != T(1)) MultXV(x,m0.diag());
      }
      inline void AssignToD(const DiagMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.size() == size());
	m2.Div(m1,m0);
	if (x != T(1)) MultXV(x,m0.diag());
      }
    protected:
      const T x;
      const GenDiagMatrix<T1>& m1;
      const GenDiagMatrix<T2>& m2;
  };

  template <class T, class T1, class T2> class RQuotDD : 
    public QuotDD<T,T1,T2>
  {
    public:
      inline RQuotDD(const T _x, const GenDiagMatrix<T1>& _m1,
	  const GenDiagMatrix<T2>& _m2) : QuotDD<T,T1,T2>(_x,_m1,_m2) {}
  };

  template <class T> inline const DiagMatrixView<T>& operator/=(
      const DiagMatrixView<T>& m1, const GenDiagMatrix<T>& m2)
  { m2.DivEq(m1); return m1; }

  template <class T> inline const DiagMatrixView<CT>& operator/=(
      const DiagMatrixView<CT>& m1, const GenDiagMatrix<T>& m2)
  { m2.DivEq(m1); return m1; }

  template <class T> inline const DiagMatrixView<T>& operator%=(
      const DiagMatrixView<T>& m1, const GenDiagMatrix<T>& m2)
  { m2.DivEq(m1); return m1; }

  template <class T> inline const DiagMatrixView<CT>& operator%=(
      const DiagMatrixView<CT>& m1, const GenDiagMatrix<T>& m2)
  { m2.DivEq(m1); return m1; }

#define QUOTMM QuotDD
#define RQUOTMM RQuotDD
#define GENMATRIX1 GenDiagMatrix
#define GENMATRIX2 GenDiagMatrix
#define PRODXM1 ProdXD
#define PRODXM2 ProdXD
#include "TMV_AuxQuotMM.h"
  // Defines things like m/m, m%m, (x*m)/m, m/(x*m), (x/m)*m, etc.
#include "TMV_AuxQuotMMa.h"
  // Defines things like -(x*m/m), (x*m/m)*x, etc.
#undef QUOTMM
#undef RQUOTMM
#undef GENMATRIX1
#undef GENMATRI2
#undef PRODXM1
#undef PRODXM2



  // Now use some helper files for the combinations which 
  // return a MatrixComposite or VectorComposite
 
  // First the VectorComposites:
  
#define GENMATRIX GenDiagMatrix
#define PRODXM ProdXD
#define PRODMV ProdDV
#define QUOTVM QuotVD
#define RQUOTVM RQuotVD
#define CONSTMATRIXVIEW ConstDiagMatrixView
 
#include "TMV_AuxVecComposite.h"

#undef PRODMV
#undef QUOTVM
#undef RQUOTVM
#undef CONSTMATRIXVIEW
#undef PRODXM
#undef GENMATRIX

  // Next all the (M) op (D) possibilities

#define GENMATRIX1 GenMatrix
#define GENMATRIX2 GenDiagMatrix
#define PRODXM1 ProdXM
#define PRODXM2 ProdXD
#define SUMMM SumMD
#define PRODMM ProdMD
#define QUOTMM QuotMD
#define RQUOTMM RQuotMD
#define TQUOTMM TransientQuotMD
#define TRQUOTMM TransientRQuotMD

#include "TMV_AuxMatComposite1.h"
#include "TMV_AuxSumMMb.h"

#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef SUMMM
#undef PRODMM
#undef QUOTXM
#undef QUOTMM
#undef RQUOTMM
#undef TQUOTMM
#undef TRQUOTMM

  // Next (D) op (M)
  
#define GENMATRIX1 GenDiagMatrix
#define GENMATRIX2 GenMatrix
#define PRODXM1 ProdXD
#define PRODXM2 ProdXM
#define SUMMMa SumMD
#define SUMMM SumDM
#define PRODMM ProdDM
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

} // namespace tmv

#undef CT
#undef CCT
#undef VCT

#ifdef TMV_TriMatrixArith_H
#include "TMV_DiagTriArith.h"
#endif

#ifdef TMV_BandMatrixArith_H
#include "TMV_DiagBandArith.h"
#endif

#ifdef TMV_SymMatrixArith_H
#include "TMV_DiagSymArith.h"
#endif

#ifdef TMV_SymBandMatrixArith_H
#include "TMV_DiagSymBandArith.h"
#endif

#endif
