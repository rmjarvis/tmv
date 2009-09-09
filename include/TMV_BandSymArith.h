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


#ifndef TMV_BandSymArith_H
#define TMV_BandSymArith_H

#define CT std::complex<T>

namespace tmv {

  //
  // BandMatrix + SymMatrix
  //

  template <class T, class T1, class T2> class SumBS : 
    public MatrixComposite<T> 
  {
    public:
      inline SumBS(T _x1, const GenBandMatrix<T1>& _m1, 
	  T _x2, const GenSymMatrix<T2>& _m2) :
	x1(_x1),m1(_m1),x2(_x2),m2(_m2)
      { 
	TMVAssert(m1.colsize() == m2.size()); 
	TMVAssert(m1.rowsize() == m2.size()); 
      }
      inline size_t colsize() const { return m2.size(); }
      inline size_t rowsize() const { return m2.size(); }
      inline StorageType stor() const { return BaseStorOf(m2); }
      inline T GetX1() const { return x1; }
      inline const GenBandMatrix<T1>& GetM1() const { return m1; }
      inline T GetX2() const { return x2; }
      inline const GenSymMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToM(const MatrixView<RealType(T)>& m0) const
      { 
	TMVAssert(IsReal(T()));
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	if (SameStorage(m0,m1)) {
	  BandMatrix<T1> m1x = m1;
	  MultXM(x2,m0=m2);
	  AddMM(x1,m1x,m0);
	} else {
	  MultXM(x2,m0=m2);
	  AddMM(x1,m1,m0);
	}
      }
      inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
      { 
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	if (SameStorage(m0,m1)) {
	  BandMatrix<T1> m1x = m1;
	  MultXM(x2,m0=m2);
	  AddMM(x1,m1x,m0);
	} else {
	  MultXM(x2,m0=m2);
	  AddMM(x1,m1,m0);
	}
      }
    private:
      const T x1;
      const GenBandMatrix<T1>& m1;
      const T x2;
      const GenSymMatrix<T2>& m2;
  };

#define SUMMM SumBS
#define GENMATRIX1 GenBandMatrix
#define GENMATRIX2 GenSymMatrix
#define PRODXM1 ProdXB
#define PRODXM2 ProdXS
#include "TMV_AuxSumMM.h"
#include "TMV_AuxSumMMb.h"
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

  //
  // BandMatrix * SymMatrix
  //

  template <class T, class T1, class T2> class ProdBS : 
    public MatrixComposite<T>
  {
    public:
      inline ProdBS(T _x, const GenBandMatrix<T1>& _m1,
	  const GenSymMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert(m1.rowsize() == m2.colsize()); }
      inline size_t colsize() const { return m1.colsize(); }
      inline size_t rowsize() const { return m2.rowsize(); }
      inline StorageType stor() const { return BaseStorOf(m2); }
      inline T GetX() const { return x; }
      inline const GenBandMatrix<T1>& GetM1() const { return m1; }
      inline const GenSymMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToM(const MatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	Matrix<T> xm2 = m2;
	MultMM<false>(x,m1,xm2,m0);
      }
      inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	Matrix<T> xm2 = m2;
	MultMM<false>(x,m1,xm2,m0);
      }
    protected:
      T x;
      const GenBandMatrix<T1>& m1;
      const GenSymMatrix<T2>& m2;
  };

  template <class T, class T1, class T2> class ProdSB : 
    public MatrixComposite<T>
  {
    public:
      inline ProdSB(T _x, const GenSymMatrix<T1>& _m1,
	  const GenBandMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert(m1.rowsize() == m2.colsize()); }
      inline size_t colsize() const { return m1.colsize(); }
      inline size_t rowsize() const { return m2.rowsize(); }
      inline StorageType stor() const { return BaseStorOf(m1); }
      inline T GetX() const { return x; }
      inline const GenSymMatrix<T1>& GetM1() const { return m1; }
      inline const GenBandMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToM(const MatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	Matrix<T> xm1 = m1;
	MultMM<false>(x,xm1,m2,m0);
      }
      inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	Matrix<T> xm1 = m1;
	MultMM<false>(x,xm1,m2,m0);
      }
    protected:
      T x;
      const GenSymMatrix<T1>& m1;
      const GenBandMatrix<T2>& m2;
  };

#define PRODMM ProdSB
#define GENMATRIX1 GenSymMatrix
#define GENMATRIX2 GenBandMatrix
#define PRODXM1 ProdXS
#define PRODXM2 ProdXB
#include "TMV_AuxProdMM.h"
#include "TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
  
#define PRODMM ProdBS
#define GENMATRIX1 GenBandMatrix
#define GENMATRIX2 GenSymMatrix
#define PRODXM1 ProdXB
#define PRODXM2 ProdXS
#include "TMV_AuxProdMM.h"
#include "TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
  
  //
  // SymMatrix / % BandMatrix
  //

#define GENMATRIX1 GenBandMatrix
#define GENMATRIX2 GenSymMatrix
#define PRODXM1 ProdXB
#define PRODXM2 ProdXS
#define QUOTXM QuotXS
#define TQUOTMM TransientQuotMS
#define TRQUOTMM TransientRQuotMS
#include "TMV_AuxTQuotMM.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTXM
#undef TQUOTMM
#undef TRQUOTMM

#define GENMATRIX1 GenSymMatrix
#define GENMATRIX2 GenBandMatrix
#define PRODXM1 ProdXS
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
