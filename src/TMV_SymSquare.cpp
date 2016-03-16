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



#include "TMV_SymSquare.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_MatrixArith.h"

namespace tmv {

    template <bool herm, class T> 
    void SymSquare(MatrixView<T> A)
    {
        const ptrdiff_t N = A.colsize();
        if (N == 1) {
            const T A00 = *A.ptr();
#ifdef TMVFLDEBUG
            TMVAssert(A.ptr() >= A._first);
            TMVAssert(A.ptr() < A._last);
#endif
            if (herm)
                *A.ptr() = TMV_NORM(TMV_REAL(A00));
            else 
                *A.ptr() = TMV_SQR(A00);
        } else {
            const ptrdiff_t K = N/2;
            MatrixView<T> A00 = A.subMatrix(0,K,0,K);
            MatrixView<T> A10 = A.subMatrix(K,N,0,K);
            MatrixView<T> A01 = A.subMatrix(0,K,K,N);
            MatrixView<T> A11 = A.subMatrix(K,N,K,N);
            MatrixView<T> A10t = herm ? A10.adjoint() : A10.transpose();

            // [ A00 A10t ] [ A00 A10t ] 
            // [ A10 A11  ] [ A10 A11  ]
            // = [ A00^2 + A10t A10    A00 A10t + A10t A11 ]
            //   [ A10 A00 + A11 A10   A10 A10t + A11^2    ]

            // A10 stores the actual data for A10
            // We can therefore write to A01 as a temp matrix.

            A01 = A00 * A10t;
            A01 += A10t * A11;

            SymSquare<herm>(A00);
            A00 += A10t*A10;
            SymSquare<herm>(A11);
            A11 += A10*A10t;

            A10t = A01;
        }
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_SymSquare.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


