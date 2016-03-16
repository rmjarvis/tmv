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

#include "TMV_SymSVDiv.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_DiagMatrixArith.h"
#include "tmv/TMV_SymMatrixArithFunc.h"
#include <ostream>

namespace tmv {


    template <class T, class T1> 
    void HermSV_Inverse(
        const GenMatrix<T1>& U, const GenDiagMatrix<TMV_RealType(T1)>& SS,
        ptrdiff_t kmax, SymMatrixView<T> sinv)
    {
        TMVAssert(sinv.isherm());
        Matrix<T,ColMajor> SinvUt = U.adjoint().rowRange(0,kmax) /
            SS.subDiagMatrix(0,kmax);
        SymMultMM<false>(T(1),U.colRange(0,kmax),SinvUt,sinv);
    }

    template <class T, class T1> 
    void SymSV_Inverse(
        const GenMatrix<T1>& U, const GenDiagMatrix<TMV_RealType(T1)>& SS, 
        const GenMatrix<T1>& Vt, ptrdiff_t kmax, SymMatrixView<T> sinv)
    {
        // A = U S Vt
        // A^-1 = V S^-1 Ut
        Matrix<T,ColMajor> SinvUt = U.adjoint().rowRange(0,kmax) /
            SS.subDiagMatrix(0,kmax);
        SymMultMM<false>(T(1),Vt.adjoint().colRange(0,kmax),SinvUt,sinv);
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_SymSVInverse.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


