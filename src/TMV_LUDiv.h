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

#ifndef TMV_LUDiv_H
#define TMV_LUDiv_H

#include "tmv/TMV_BaseMatrix.h"

namespace tmv {

    template <typename T, typename T1> 
    void LU_LDivEq(
        const GenMatrix<T1>& LUx, const ptrdiff_t* P, MatrixView<T> m);

    template <typename T, typename T1> 
    void LU_RDivEq(
        const GenMatrix<T1>& LUx, const ptrdiff_t* P, MatrixView<T> m);

    template <typename T, typename T1> 
    void LU_Inverse(
        const GenMatrix<T1>& LUx, const ptrdiff_t* P, MatrixView<T> m);

    // Specialize disallowed complex combinations:
#define CT std::complex<T>
    template <typename T> 
    inline void LU_LDivEq(
        const GenMatrix<CT>& , const ptrdiff_t* , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <typename T> 
    inline void LU_RDivEq(
        const GenMatrix<CT>& , const ptrdiff_t* , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <typename T> 
    inline void LU_Inverse(
        const GenMatrix<CT>& , const ptrdiff_t* , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
#undef CT

}

#endif
