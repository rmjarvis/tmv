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


// This file sets up the Composite QuotXM classes for Matrix types which
// returns a Matrix for the operations x/m or x%m.
// (Some do not: DiagMatrix for example, returns another DiagMatix.)
//
// Need to define the following with #define statements.
// (The given definition is for a Band Matrix.  Modify as 
// appropriate for the various other matrices.)
//
// #define GENMATRIX2 GenBandMatrix
//
// #define QUOTXM QuotXB


//
// Scalar / Matrix
//

template <class T, class Tm> class QUOTXM : 
  public MatrixComposite<T> 
{
  public:
    inline QUOTXM(const T _x, const GENMATRIX2<Tm>& _m) : x(_x), m(_m) {}
    inline size_t colsize() const { return m.rowsize(); }
    inline size_t rowsize() const { return m.colsize(); }
    inline StorageType stor() const { return ColMajor; }
    inline T GetX() const { return x; }
    inline const GENMATRIX2<Tm>& GetM() const { return m; }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    {
      TMVAssert(IsReal(T()));
      TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
      m.Inverse(m0);
      MultXM(x,m0);
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    {
      TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
      m.Inverse(m0);
      MultXM(x,m0);
    }
  private:
    const T x;
    const GENMATRIX2<Tm>& m;
};

#include "TMV_AuxQuotXM.h"
