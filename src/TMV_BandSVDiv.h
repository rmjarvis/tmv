///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2014                                                 //
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

#ifndef TMV_BandSVDiv_H
#define TMV_BandSVDiv_H

#include "tmv/TMV_BaseBandMatrix.h"
#include "TMV_SVDiv.h"

namespace tmv {

    template <class T> 
    void SV_Decompose(
        const GenBandMatrix<T>& A,
        MatrixView<T> U, DiagMatrixView<TMV_RealType(T)> S, MatrixView<T> Vt, 
        TMV_RealType(T)& logdet, T& signdet);

}

#endif
