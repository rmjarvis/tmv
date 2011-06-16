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
#include "tmv/TMV_QRInverse.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_CopyV.h"
#include "tmv/TMV_DivVU.h"
#include "tmv/TMV_DivMU.h"
#include "tmv/TMV_PermuteM.h"
#include "tmv/TMV_MultUL.h"
#include "tmv/TMV_MultPM.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_ScaleM.h"
#include "tmv/TMV_ProdMM.h"
#include "tmv/TMV_QuotMM.h"

namespace tmv {

    //
    // Inverse
    //
    
    template <class M1, class V1, class M2>
    static void DoQR_Inverse(
        const M1& QR, const V1& beta, const Permutation* P, int N1, M2& m2)
    {
        if (beta.step() == 1) {
            if (QR.iscm()) {
                InlineQR_Inverse(
                    QR.cmView(),beta.unitView(),P,N1,m2);
            } else if (QR.isrm()) {
                InlineQR_Inverse(
                    QR.rmView(),beta.unitView(),P,N1,m2);
            } else {
                InlineQR_Inverse(
                    QR.copy().view(),beta.unitView(),P,N1,m2);
            }
        } else {
            DoQR_Inverse(QR,beta.copy().xView(),P,N1,m2);
        }
    }

    template <class T1, int C1, class RT1, class T2>
    void InstQR_Inverse(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1, MatrixView<T2> m2)
    {
        if (m2.iscm()) {
            MatrixView<T2,ColMajor> m2cm = m2;
            DoQR_Inverse(QR,beta,P,N1,m2cm);
        } else if (m2.isrm()) {
            MatrixView<T2,RowMajor> m2rm = m2;
            DoQR_Inverse(QR,beta,P,N1,m2rm);
        } else {
            Matrix<T2,ColMajor|NoDivider> m2c(m2);
            MatrixView<T2,ColMajor> m2cm = m2c.view();
            DoQR_Inverse(QR,beta,P,N1,m2cm);
            InstCopy(m2cm.constView().xView(),m2);
        }
    }

    //
    // InverseATA
    //

    template <class M1, class V1, class M2>
    static void DoQR_InverseATA(
        const M1& QR, const V1& beta, const Permutation* P, int N1, M2& m2)
    {
        if (beta.step() == 1) {
            if (QR.iscm()) {
                InlineQR_InverseATA(
                    QR.cmView(),beta.unitView(),P,N1,m2);
            } else if (QR.isrm()) {
                InlineQR_InverseATA(
                    QR.rmView(),beta.unitView(),P,N1,m2);
            } else {
                InlineQR_InverseATA(
                    QR.copy().view(),beta.unitView(),P,N1,m2);
            }
        } else {
            DoQR_InverseATA(QR,beta.copy().xView(),P,N1,m2);
        }
    }

    template <class T1, int C1, class RT1, class T2>
    void InstQR_InverseATA(
        const ConstMatrixView<T1,C1>& QR, const ConstVectorView<RT1>& beta,
        const Permutation* P, int N1, MatrixView<T2> m2)
    {
        if (m2.iscm()) {
            MatrixView<T2,ColMajor> m2cm = m2;
            DoQR_InverseATA(QR,beta,P,N1,m2cm);
        } else if (m2.isrm()) {
            MatrixView<T2,RowMajor> m2rm = m2;
            DoQR_InverseATA(QR,beta,P,N1,m2rm);
        } else {
            Matrix<T2,ColMajor|NoDivider> m2c(m2);
            MatrixView<T2,ColMajor> m2cm = m2c.view();
            DoQR_InverseATA(QR,beta,P,N1,m2cm);
            InstCopy(m2cm.constView().xView(),m2);
        }
    }




#define InstFile "TMV_QRInverse.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


