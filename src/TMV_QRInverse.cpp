///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2016                                                 //
// All rights reserved                                                       //
//                                                                           //
// The project is hosted at https://code.google.com/p/tmv-cpp/               //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis17 [at] gmail.                                                 //
//                                                                           //
// This software is licensed under a FreeBSD license.  The file              //
// TMV_LICENSE should have bee included with this distribution.              //
// It not, you can get a copy from https://code.google.com/p/tmv-cpp/.       //
//                                                                           //
// Essentially, you can use this software however you want provided that     //
// you include the TMV_LICENSE file in any distribution that uses it.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TMV_QRDiv.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_TriMatrix.h"

namespace tmv {

    template <class T, class T1> 
    void QR_Inverse(
        const GenMatrix<T1>& QRx, const GenVector<T1>& beta, const ptrdiff_t* P,
        MatrixView<T> minv, ptrdiff_t N1)
    {
        // minv = Pt R^-1 Qt
        TMVAssert(minv.colsize() == QRx.rowsize());
        TMVAssert(minv.rowsize() == QRx.colsize());

        minv.setZero();
        UpperTriMatrixView<T> R = minv.colRange(0,N1).upperTri();
        R = QRx.upperTri().subTriMatrix(0,N1);
        R.invertSelf();
        Q_RDivEq(QRx,beta,minv);
        if (P) minv.reversePermuteRows(P);
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_QRInverse.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


