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


// This file sets up the Composite classes for all operations with a 
// (sparse) matrix that returns a Matrix and is of the form:
// (S) op (M)
// (ie. the second operand is the normal Matrix.)
//
// Need to define the following with #define statements.
// (The given definition is for a Band Matrix.  Modify as 
// appropriate for the various other matrices.)
//
// #define GENMATRIX1 GenBandMatrix
//
// #define PRODMM ProdBM


//
// Matrix * Matrix
//

template <class T, class T1, class T2> class PRODMM : 
  public MatrixComposite<T>
{
  public:
    inline PRODMM(const T _x, const GENMATRIX1<T1>& _m1,
	const GenMatrix<T2>& _m2) :
      x(_x), m1(_m1), m2(_m2)
    { TMVAssert(m1.rowsize() == m2.colsize()) ; }
    inline size_t colsize() const { return m1.colsize(); }
    inline size_t rowsize() const { return m2.rowsize(); }
    inline StorageType stor() const { return BaseStorOf(m2); }
    inline T GetX() const { return x; }
    inline const GENMATRIX1<T1>& GetM1() const { return m1; }
    inline const GenMatrix<T2>& GetM2() const { return m2; }
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
    const GENMATRIX1<T1>& m1;
    const GenMatrix<T2>& m2;
};

template <class T, class T2, class T3> inline const MatrixView<T>& operator+=(
    const MatrixView<T>& m, const PRODMM<T,T2,T3>& pmm)
{
  TMVAssert(m.colsize() == pmm.colsize());
  TMVAssert(m.rowsize() == pmm.rowsize());
  MultMM<true>(pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m);
  return m;
}

template <class T> inline const MatrixView<CT>& operator+=(
    const MatrixView<CT>& m, const PRODMM<T,T,T>& pmm)
{
  TMVAssert(m.colsize() == pmm.colsize());
  TMVAssert(m.rowsize() == pmm.rowsize());
  MultMM<true>(CT(pmm.GetX()),pmm.GetM1(),pmm.GetM2(),m);
  return m;
}

template <class T, class T2, class T3> inline const MatrixView<T>& operator-=(
    const MatrixView<T>& m, const PRODMM<T,T2,T3>& pmm)
{
  TMVAssert(m.colsize() == pmm.colsize());
  TMVAssert(m.rowsize() == pmm.rowsize());
  MultMM<true>(-pmm.GetX(),pmm.GetM1(),pmm.GetM2(),m);
  return m;
}

template <class T> inline const MatrixView<CT>& operator-=(
    const MatrixView<CT>& m, const PRODMM<T,T,T>& pmm)
{
  TMVAssert(m.colsize() == pmm.colsize());
  TMVAssert(m.rowsize() == pmm.rowsize());
  MultMM<true>(CT(-pmm.GetX()),pmm.GetM1(),pmm.GetM2(),m);
  return m;
}

#include "TMV_AuxProdMM.h"
#include "TMV_AuxProdMMa.h"


