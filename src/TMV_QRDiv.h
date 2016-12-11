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

#ifndef TMV_QRDiv_H
#define TMV_QRDiv_H

#include "tmv/TMV_BaseMatrix.h"
#include "tmv/TMV_BaseTriMatrix.h"

namespace tmv {

    // Decompose A (input as QRx) into Q R.
    // On output, Q is stored in the lower triangle part of QRx as
    // Householder vectors, and the beta vector.
    // R is the upper triangle part of QRx
    template <typename T> 
    void QR_Decompose(
        MatrixView<T> QRx, VectorView<T> beta, T& signdet);

    template <typename T> 
    void QR_Decompose(
        MatrixView<T> Q, UpperTriMatrixView<T> R, T& signdet);

    template <typename T> 
    void GetQFromQR(MatrixView<T> Q, const GenVector<T>& beta);

    template <typename T, typename T1> 
    void Q_LDivEq(
        const GenMatrix<T1>& Q, const GenVector<T1>& beta, 
        MatrixView<T> m);
    template <typename T, typename T1> 
    void Q_RDivEq(
        const GenMatrix<T1>& Q, const GenVector<T1>& beta, 
        MatrixView<T> m);


    template <typename T, typename T1, typename T2> 
    void QR_LDiv(
        const GenMatrix<T1>& QR, const GenVector<T1>& beta, const ptrdiff_t* P,
        const GenMatrix<T2>& m, MatrixView<T> x, ptrdiff_t N1);
    template <typename T, typename T1> 
    void QR_LDivEq(
        const GenMatrix<T1>& QR, const GenVector<T1>& beta, const ptrdiff_t* P,
        MatrixView<T> m, ptrdiff_t N1);
    template <typename T, typename T1, typename T2> 
    void QR_RDiv(
        const GenMatrix<T1>& QR, const GenVector<T1>& beta, const ptrdiff_t* P,
        const GenMatrix<T2>& m, MatrixView<T> x, ptrdiff_t N1);
    template <typename T, typename T1> 
    void QR_RDivEq(
        const GenMatrix<T1>& QR, const GenVector<T1>& beta, const ptrdiff_t* P,
        MatrixView<T> m, ptrdiff_t N1);
    template <typename T, typename T1> 
    void QR_Inverse(
        const GenMatrix<T1>& QRx, const GenVector<T1>& beta, const ptrdiff_t* P,
        MatrixView<T> minv, ptrdiff_t N1);

    // Specialize disallowed complex combinations:
#define CT std::complex<T>
    template <typename T>
    inline void Q_LDivEq(
        const GenMatrix<CT>& , const GenVector<CT>& , 
        MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T>
    inline void Q_RDivEq(
        const GenMatrix<CT>& , const GenVector<CT>& , 
        MatrixView<T> )
    { TMVAssert(TMV_FALSE); }


    template <typename T>
    inline void QR_LDiv(
        const GenMatrix<CT>& , const GenVector<CT>& , const ptrdiff_t* ,
        const GenMatrix<CT>& , MatrixView<T> , ptrdiff_t )
    { TMVAssert(TMV_FALSE); }
    template <typename T>
    inline void QR_LDiv(
        const GenMatrix<CT>& , const GenVector<CT>& , const ptrdiff_t* ,
        const GenMatrix<T>& , MatrixView<T> , ptrdiff_t )
    { TMVAssert(TMV_FALSE); }
    template <typename T>
    inline void QR_LDiv(
        const GenMatrix<T>& , const GenVector<T>& , const ptrdiff_t* ,
        const GenMatrix<CT>& , MatrixView<T> , ptrdiff_t )
    { TMVAssert(TMV_FALSE); }
    template <typename T>
    inline void QR_LDivEq(
        const GenMatrix<CT>& , const GenVector<CT>& , const ptrdiff_t* ,
        MatrixView<T> , ptrdiff_t )
    { TMVAssert(TMV_FALSE); }
    template <typename T>
    inline void QR_RDiv(
        const GenMatrix<CT>& , const GenVector<CT>& , const ptrdiff_t* ,
        const GenMatrix<CT>& , MatrixView<T> , ptrdiff_t )
    { TMVAssert(TMV_FALSE); }
    template <typename T>
    inline void QR_RDiv(
        const GenMatrix<CT>& , const GenVector<CT>& , const ptrdiff_t* ,
        const GenMatrix<T>& , MatrixView<T> , ptrdiff_t )
    { TMVAssert(TMV_FALSE); }
    template <typename T>
    inline void QR_RDiv(
        const GenMatrix<T>& , const GenVector<T>& , const ptrdiff_t* ,
        const GenMatrix<CT>& , MatrixView<T> , ptrdiff_t )
    { TMVAssert(TMV_FALSE); }
    template <typename T>
    inline void QR_RDivEq(
        const GenMatrix<CT>& , const GenVector<CT>& , const ptrdiff_t* ,
        MatrixView<T> , ptrdiff_t )
    { TMVAssert(TMV_FALSE); }
    template <typename T>
    inline void QR_Inverse(
        const GenMatrix<CT>& , const GenVector<CT>& , const ptrdiff_t* ,
        MatrixView<T> , ptrdiff_t )
    { TMVAssert(TMV_FALSE); }
#undef CT


}

#endif
