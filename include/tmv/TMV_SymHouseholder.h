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


#ifndef TMV_SymHouseholder_H
#define TMV_SymHouseholder_H

#include "tmv/TMV_BaseSymMatrix.h"

namespace tmv {

    template <typename T1, typename T2>
    void HouseholderLRMult(
        const GenVector<T1>& v, T1 beta, SymMatrixView<T2> M);
    // The input vector, v, is taken to be the vector for a
    // Householder matrix, H.  This routine takes M <- H M Ht
    // if M is Hermitian or H M HT if M is symmetric.

} // namespace tmv

#endif
