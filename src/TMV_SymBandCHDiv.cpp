///////////////////////////////////////////////////////////////////////////////
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



#include "TMV_SymBandCHDiv.h"
#include "tmv/TMV_SymBandMatrix.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_DiagMatrixArith.h"
#include "TMV_BandLUDiv.h"
#include "tmv/TMV_BandMatrixArith.h"
#include <ostream>

namespace tmv {

    template <class T, class T1> 
    void CH_LDivEq(const GenSymBandMatrix<T1>& LL, MatrixView<T> m)
    {
        TMVAssert(LL.size() == m.colsize());
        // m = (LLt)^-1 m
        //   = Lt^-1 L^-1 m
        m /= LL.lowerBand();
        m /= LL.upperBand();
    }

    template <class T, class T1> 
    void LDL_LDivEq(const GenSymBandMatrix<T1>& LL, MatrixView<T> m)
    {
        TMVAssert(LL.size() == m.colsize());
        TMVAssert(LL.nlo() == 1);
        // m = (LDLt)^-1 m
        //   = Lt^-1 D^-1 L^-1 m

        // m /= LL.lowerBand(UnitDiag)
        TriLDivEq(LL.lowerBand(),m,UnitDiag);
        m /= DiagMatrixViewOf(LL.diag().realPart());
        // m /= LL.upperBand(UnitDiag)
        TriLDivEq(LL.upperBand(),m,UnitDiag);
    }

    template <class T, class T1> 
    void CH_RDivEq(const GenSymBandMatrix<T1>& LL, MatrixView<T> m)
    {
        TMVAssert(LL.size() == m.rowsize());
        // m = m (LLt)^-1
        //   = m Lt^-1 L^-1
        // mt = Lt^-1 L^-1 mt
        m.adjoint() /= LL.lowerBand();
        m.adjoint() /= LL.upperBand();
    }

    template <class T, class T1> 
    void LDL_RDivEq(const GenSymBandMatrix<T1>& LL, MatrixView<T> m)
    {
        TMVAssert(LL.size() == m.rowsize());
        TMVAssert(LL.nlo() == 1);
        // m = m (LDLt)^-1
        //   = m Lt^-1 D^-1 L^-1
        // mt = Lt^-1 D^-1 L^-1 mt

        // mt /= LL.lowerBand(UnitDiag)
        TriLDivEq(LL.lowerBand(),m.adjoint(),UnitDiag);
        m %= DiagMatrixViewOf(LL.diag().realPart());
        // mt /= LL.upperBand(UnitDiag)
        TriLDivEq(LL.upperBand(),m.adjoint(),UnitDiag);
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_SymBandCHDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


