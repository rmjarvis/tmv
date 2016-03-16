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


#include "TMV_SymCHDiv.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_TriMatrixArith.h"
#include <ostream>

namespace tmv {

    //
    // LDivEq
    //

    template <class T, class T1> 
    void CH_LDivEq(const GenSymMatrix<T1>& LL, MatrixView<T> m)
    {
        TMVAssert(LL.size() == m.colsize());
        // m = (LLt)^-1 m
        //   = Lt^-1 L^-1 m
        m /= LL.lowerTri();
        m /= LL.upperTri();
    }

    //
    // RDivEq Matrix
    //

    template <class T, class T1> 
    void CH_RDivEq(const GenSymMatrix<T1>& LL, MatrixView<T> m)
    {
        TMVAssert(LL.size() == m.rowsize());
        // m = m (LLt)^-1 
        //   = m Lt^-1 L^-1
        m %= LL.upperTri();
        m %= LL.lowerTri();
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_SymCHDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


