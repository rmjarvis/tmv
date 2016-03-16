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



#include "tmv/TMV_DiagMatrixArithFunc.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_MatrixArith.h"

namespace tmv {

    template <class T, class Ta, class Tb> 
    void AddMM(
        const T alpha, const GenDiagMatrix<Ta>& A,
        const T beta, const GenMatrix<Tb>& B, MatrixView<T> C)
        // C = alpha*A + beta*B
    {
        TMVAssert(A.size() == B.colsize());
        TMVAssert(A.size() == B.rowsize());
        TMVAssert(A.size() == C.colsize());
        TMVAssert(A.size() == C.rowsize());

        if (A.size() > 0) {
            if (SameStorage(A.diag(),C)) {
                DiagMatrix<Ta> tempA = A;
                C = beta * B;
                AddMM(alpha,tempA,C);
            } else {
                C = beta * B;
                AddMM(alpha,A,C);
            }
        }
    }

#define InstFile "TMV_AddDM.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


