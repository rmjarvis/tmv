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

//#define PRINTALGO_QR

#include "TMV_Blas.h"
#include "tmv/TMV_QRDecompose.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_CopyV.h"
#include "tmv/TMV_SwapV.h"
#include "tmv/TMV_MinMax.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_Rank1VVM.h"
#include "tmv/TMV_MultMV.h"
#include "tmv/TMV_NormM.h"
#include "tmv/TMV_MultMM.h"

namespace tmv {

#ifdef LAP
    template <class T, int S, class RT> 
    static inline void LapQR_Decompose(
        MatrixView<T,S> A, VectorView<RT,Unit> beta)
    { InlineQR_Decompose(A,beta); }
#endif // LAP

    template <class T, class RT> 
    void InstQR_Decompose(MatrixView<T> A, VectorView<RT> beta)
    {
        if (A.rowsize() > 0) {
            if (beta.step() == 1) {
                VectorView<RT,Unit> beta1 = beta;
                if (A.iscm()) {
                    MatrixView<T,ColMajor> Acm = A;
#ifdef LAP
                    LapQR_Decompose(Acm,beta1);
#else
                    InlineQR_Decompose(Acm,beta1);
#endif
                } else if (A.isrm()) {
                    MatrixView<T,RowMajor> Arm = A;
#ifdef LAP
                    LapQR_Decompose(Arm,beta1);
#else
                    InlineQR_Decompose(Arm,beta1);
#endif
                } else {
                    Matrix<T,ColMajor|NoDivider> Ac = A;
                    InstQR_Decompose(Ac.xView(),beta);
                    InstCopy(Ac.constView().xView(),A);
                }
            } else {
                Vector<RT> betac = beta;
                InstQR_Decompose(A,betac.xView());
                InstCopy(betac.constView().xView(),beta);
            }
        }
    }

#define InstFile "TMV_QRDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


