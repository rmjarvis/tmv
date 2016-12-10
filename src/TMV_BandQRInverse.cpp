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


#include "TMV_BandQRDiv.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "TMV_BandLUDiv.h"
#include "TMV_BandQRDiv.h"

namespace tmv {

    template <class T, class T1> 
    void QR_Inverse(
        const GenBandMatrix<T1>& QRx, const GenVector<T1>& beta,
        const MatrixView<T>& minv)
    {
        // minv = R^-1 Qt
        TMVAssert(minv.colsize() == QRx.rowsize());
        TMVAssert(minv.rowsize() == QRx.colsize());

        const ptrdiff_t N = QRx.rowsize();

        minv.setZero();
        UpperTriMatrixView<T> R = minv.colRange(0,N).upperTri();
        R = QRx.diagRange(0,QRx.nhi()+1);
        TriInverse(R,QRx.nhi());
        Q_RDivEq(QRx,beta,minv);
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_BandQRInverse.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


