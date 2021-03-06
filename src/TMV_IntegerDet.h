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

#ifndef TMV_IntegerDet_H
#define TMV_IntegerDet_H

namespace tmv {

    template <typename T> 
    T IntegerDet(const GenMatrix<T>& A);
    template <typename T> 
    T IntegerDet(const GenBandMatrix<T>& A);
    template <typename T> 
    T IntegerDet(const GenSymMatrix<T>& A);
    template <typename T> 
    T IntegerDet(const GenSymBandMatrix<T>& A);

} // namespace tmv

#endif
