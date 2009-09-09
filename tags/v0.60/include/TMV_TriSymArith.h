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


#ifndef TMV_TriSymArith_H
#define TMV_TriSymArith_H

#define CT std::complex<T>

namespace tmv {

  //
  // TriMatrix + SymMatrix
  //

  template <class T, class T1, class T2> class SumUS : 
    public MatrixComposite<T> 
  {
    public:
      inline SumUS(T _x1, const GenUpperTriMatrix<T1>& _m1, 
	  T _x2, const GenSymMatrix<T2>& _m2) :
	x1(_x1),m1(_m1),x2(_x2),m2(_m2)
      { TMVAssert(m1.size() == m2.size()); }
      inline size_t colsize() const { return m1.size(); }
      inline size_t rowsize() const { return m1.size(); }
      inline StorageType stor() const { return BaseStorOf(m2); }
      inline T GetX1() const { return x1; }
      inline const GenUpperTriMatrix<T1>& GetM1() const { return m1; }
      inline T GetX2() const { return x2; }
      inline const GenSymMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToM(const MatrixView<RealType(T)>& m0) const
      { 
	TMVAssert(IsReal(T()));
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	if (SameStorage(m0,m1)) {
	  UpperTriMatrix<T1> m1x = m1;
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
	  UpperTriMatrix<T1> m1x = m1;
	  MultXM(x2,m0=m2);
	  AddMM(x1,m1x,m0);
	} else {
	  MultXM(x2,m0=m2);
	  AddMM(x1,m1,m0);
	}
      }
    private:
      const T x1;
      const GenUpperTriMatrix<T1>& m1;
      const T x2;
      const GenSymMatrix<T2>& m2;
  };

  template <class T, class T1, class T2> class SumLS : 
    public MatrixComposite<T> 
  {
    public:
      inline SumLS(T _x1, const GenLowerTriMatrix<T1>& _m1, 
	  T _x2, const GenSymMatrix<T2>& _m2) :
	x1(_x1),m1(_m1),x2(_x2),m2(_m2)
      { TMVAssert(m1.size() == m2.size()); }
      inline size_t colsize() const { return m1.size(); }
      inline size_t rowsize() const { return m1.size(); }
      inline StorageType stor() const { return BaseStorOf(m2); }
      inline T GetX1() const { return x1; }
      inline const GenLowerTriMatrix<T1>& GetM1() const { return m1; }
      inline T GetX2() const { return x2; }
      inline const GenSymMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToM(const MatrixView<RealType(T)>& m0) const
      { 
	TMVAssert(IsReal(T()));
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	if (SameStorage(m0,m1)) {
	  LowerTriMatrix<T1> m1x = m1;
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
	  LowerTriMatrix<T1> m1x = m1;
	  MultXM(x2,m0=m2);
	  AddMM(x1,m1x,m0);
	} else {
	  MultXM(x2,m0=m2);
	  AddMM(x1,m1,m0);
	}
      }
    private:
      const T x1;
      const GenLowerTriMatrix<T1>& m1;
      const T x2;
      const GenSymMatrix<T2>& m2;
  };

#define SUMMM SumUS
#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenSymMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXS
#include "TMV_AuxSumMM.h"
#include "TMV_AuxSumMMb.h"
#undef SUMMM
#undef GENMATRIX1
#undef PRODXM1
#define SUMMM SumLS
#define GENMATRIX1 GenLowerTriMatrix
#define PRODXM1 ProdXL
#include "TMV_AuxSumMM.h"
#include "TMV_AuxSumMMb.h"
#undef SUMMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2

  //
  // TriMatrix * SymMatrix
  //

  template <class T, class T1, class T2> class ProdUS : 
    public MatrixComposite<T>
  {
    public:
      inline ProdUS(T _x, const GenUpperTriMatrix<T1>& _m1,
	  const GenSymMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert(m1.size() == m2.size()); }
      inline size_t colsize() const { return m2.size(); }
      inline size_t rowsize() const { return m2.size(); }
      inline StorageType stor() const { return BaseStorOf(m2); }
      inline T GetX() const { return x; }
      inline const GenUpperTriMatrix<T1>& GetM1() const { return m1; }
      inline const GenSymMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToM(const MatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	if (SameStorage(m0,m1)) {
	  UpperTriMatrix<T1> m1x = m1;
	  MultXM(x,m0=m2);
	  MultMM<false>(T(1),m1x,m0,m0);
	} else {
	  MultXM(x,m0=m2);
	  MultMM<false>(T(1),m1,m0,m0);
	}
      }
      inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	if (SameStorage(m0,m1)) {
	  UpperTriMatrix<T1> m1x = m1;
	  MultXM(x,m0=m2);
	  MultMM<false>(T(1),m1x,m0,m0);
	} else {
	  MultXM(x,m0=m2);
	  MultMM<false>(T(1),m1,m0,m0);
	}
      }
    protected:
      T x;
      const GenUpperTriMatrix<T1>& m1;
      const GenSymMatrix<T2>& m2;
  };

  template <class T, class T1, class T2> class ProdSU : 
    public MatrixComposite<T>
  {
    public:
      inline ProdSU(T _x, const GenSymMatrix<T1>& _m1,
	  const GenUpperTriMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert(m1.size() == m2.size()); }
      inline size_t colsize() const { return m1.size(); }
      inline size_t rowsize() const { return m1.size(); }
      inline StorageType stor() const { return BaseStorOf(m1); }
      inline T GetX() const { return x; }
      inline const GenSymMatrix<T1>& GetM1() const { return m1; }
      inline const GenUpperTriMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToM(const MatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	if (SameStorage(m0,m2)) {
	  UpperTriMatrix<T2> m2x = m2;
	  MultXM(x,m0=m1);
	  MultMM<false>(T(1),m0,m2x,m0);
	} else {
	  MultXM(x,m0=m1);
	  MultMM<false>(T(1),m0,m2,m0);
	}
      }
      inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	if (SameStorage(m0,m2)) {
	  UpperTriMatrix<T2> m2x = m2;
	  MultXM(x,m0=m1);
	  MultMM<false>(T(1),m0,m2x,m0);
	} else {
	  MultXM(x,m0=m1);
	  MultMM<false>(T(1),m0,m2,m0);
	}
      }
    protected:
      T x;
      const GenSymMatrix<T1>& m1;
      const GenUpperTriMatrix<T2>& m2;
  };

  template <class T, class T1, class T2> class ProdLS : 
    public MatrixComposite<T>
  {
    public:
      inline ProdLS(T _x, const GenLowerTriMatrix<T1>& _m1,
	  const GenSymMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert(m1.size() == m2.size()); }
      inline size_t colsize() const { return m2.size(); }
      inline size_t rowsize() const { return m2.size(); }
      inline StorageType stor() const { return BaseStorOf(m2); }
      inline T GetX() const { return x; }
      inline const GenLowerTriMatrix<T1>& GetM1() const { return m1; }
      inline const GenSymMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToM(const MatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	if (SameStorage(m0,m1)) {
	  LowerTriMatrix<T1> m1x = m1;
	  MultXM(x,m0=m2);
	  MultMM<false>(T(1),m1x,m0,m0);
	} else {
	  MultXM(x,m0=m2);
	  MultMM<false>(T(1),m1,m0,m0);
	}
      }
      inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	if (SameStorage(m0,m1)) {
	  LowerTriMatrix<T1> m1x = m1;
	  MultXM(x,m0=m2);
	  MultMM<false>(T(1),m1x,m0,m0);
	} else {
	  MultXM(x,m0=m2);
	  MultMM<false>(T(1),m1,m0,m0);
	}
      }
    protected:
      T x;
      const GenLowerTriMatrix<T1>& m1;
      const GenSymMatrix<T2>& m2;
  };

  template <class T, class T1, class T2> class ProdSL : 
    public MatrixComposite<T>
  {
    public:
      inline ProdSL(T _x, const GenSymMatrix<T1>& _m1,
	  const GenLowerTriMatrix<T2>& _m2) :
	x(_x), m1(_m1), m2(_m2)
      { TMVAssert(m1.size() == m2.size()); }
      inline size_t colsize() const { return m1.size(); }
      inline size_t rowsize() const { return m1.size(); }
      inline StorageType stor() const { return BaseStorOf(m1); }
      inline T GetX() const { return x; }
      inline const GenSymMatrix<T1>& GetM1() const { return m1; }
      inline const GenLowerTriMatrix<T2>& GetM2() const { return m2; }
      inline void AssignToM(const MatrixView<RealType(T)>& m0) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	if (SameStorage(m0,m2)) {
	  LowerTriMatrix<T2> m2x = m2;
	  MultXM(x,m0=m1);
	  MultMM<false>(T(1),m0,m2x,m0);
	} else {
	  MultXM(x,m0=m1);
	  MultMM<false>(T(1),m0,m2,m0);
	}
      }
      inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.colsize() == colsize());
	TMVAssert(m0.rowsize() == rowsize());
	if (SameStorage(m0,m2)) {
	  LowerTriMatrix<T2> m2x = m2;
	  MultXM(x,m0=m1);
	  MultMM<false>(T(1),m0,m2x,m0);
	} else {
	  MultXM(x,m0=m1);
	  MultMM<false>(T(1),m0,m2,m0);
	}
      }
    protected:
      T x;
      const GenSymMatrix<T1>& m1;
      const GenLowerTriMatrix<T2>& m2;
  };

#define PRODMM ProdSU
#define GENMATRIX1 GenSymMatrix
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM1 ProdXS
#define PRODXM2 ProdXU
#include "TMV_AuxProdMM.h"
#include "TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX2
#undef PRODXM2
#define PRODMM ProdSL
#define GENMATRIX2 GenLowerTriMatrix
#define PRODXM2 ProdXL
#include "TMV_AuxProdMM.h"
#include "TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
  
#define PRODMM ProdUS
#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenSymMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXS
#include "TMV_AuxProdMM.h"
#include "TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef PRODXM1
#define PRODMM ProdLS
#define GENMATRIX1 GenLowerTriMatrix
#define PRODXM1 ProdXL
#include "TMV_AuxProdMM.h"
#include "TMV_AuxProdMMa.h"
#undef PRODMM
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
  
  //
  // SymMatrix / % TriMatrix
  //

#define GENMATRIX1 GenSymMatrix
#define GENMATRIX2 GenUpperTriMatrix
#define PRODXM1 ProdXS
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

#define GENMATRIX1 GenUpperTriMatrix
#define GENMATRIX2 GenSymMatrix
#define PRODXM1 ProdXU
#define PRODXM2 ProdXS
#define QUOTXM QuotXS
#define TQUOTMM TransientQuotMS
#define TRQUOTMM TransientRQuotMS
#include "TMV_AuxTQuotMM.h"
#undef GENMATRIX1
#undef PRODXM1
#define GENMATRIX1 GenLowerTriMatrix
#define PRODXM1 ProdXL
#include "TMV_AuxTQuotMM.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTXM
#undef TQUOTMM
#undef TRQUOTMM

}

#undef CT

#endif
