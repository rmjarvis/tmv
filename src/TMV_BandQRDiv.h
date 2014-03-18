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

#ifndef TMV_BandQRDiv_H
#define TMV_BandQRDiv_H

#include "tmv/TMV_BaseBandMatrix.h"

namespace tmv {

    template <class T> 
    void QR_Decompose(
        BandMatrixView<T> QRx, VectorView<T> Qbeta, T& signdet);

    template <class T> 
    void GetQFromBandQR(
        MatrixView<T> QRx, const GenVector<T>& Qbeta, ptrdiff_t nlo);

    template <class T, class T1> 
    void QR_LDivEq(
        const GenBandMatrix<T1>& QRx, const GenVector<T1>& Qbeta,
        MatrixView<T> m);
    template <class T, class T1> 
    void QR_RDivEq(
        const GenBandMatrix<T1>& QRx, const GenVector<T1>& Qbeta,
        MatrixView<T> m);

    template <class T, class T1, class T2> 
    void QR_LDiv(
        const GenBandMatrix<T1>& QR, const GenVector<T1>& Qbeta,
        const GenMatrix<T2>& m, MatrixView<T> x);
    template <class T, class T1, class T2> 
    void QR_RDiv(
        const GenBandMatrix<T1>& QR, const GenVector<T1>& Qbeta,
        const GenMatrix<T2>& m, MatrixView<T> x);

    template <class T, class T1> 
    void Q_LDivEq(
        const GenBandMatrix<T1>& Q, const GenVector<T1>& Qbeta,
        MatrixView<T> m);
    template <class T, class T1> 
    void Q_RDivEq(
        const GenBandMatrix<T1>& Q, const GenVector<T1>& Qbeta,
        MatrixView<T> m);

    template <class T, class T1> 
    void QR_Inverse(
        const GenBandMatrix<T1>& QRx, const GenVector<T1>& Qbeta,
        MatrixView<T> m);


    // Specialize disallowed complex combinations:
#define CT std::complex<T>

    template <class T>
    inline void QR_LDivEq(
        const GenBandMatrix<CT>& , const GenVector<CT>& ,
        MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T>
    inline void QR_RDivEq(
        const GenBandMatrix<CT>& , const GenVector<CT>& ,
        MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <class T>
    inline void QR_LDiv(
        const GenBandMatrix<CT>& , const GenVector<CT>& ,
        const GenMatrix<CT>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T>
    inline void QR_LDiv(
        const GenBandMatrix<CT>& , const GenVector<CT>& ,
        const GenMatrix<T>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T>
    inline void QR_LDiv(
        const GenBandMatrix<T>& , const GenVector<T>& ,
        const GenMatrix<CT>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T>
    inline void QR_RDiv(
        const GenBandMatrix<CT>& , const GenVector<CT>& ,
        const GenMatrix<CT>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T>
    inline void QR_RDiv(
        const GenBandMatrix<CT>& , const GenVector<CT>& ,
        const GenMatrix<T>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T>
    inline void QR_RDiv(
        const GenBandMatrix<T>& , const GenVector<T>& ,
        const GenMatrix<CT>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <class T>
    inline void Q_LDivEq(
        const GenBandMatrix<CT>& , const GenVector<CT>& ,
        MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T>
    inline void Q_RDivEq(
        const GenBandMatrix<CT>& , const GenVector<CT>& ,
        MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <class T>
    inline void QR_Inverse(
        const GenBandMatrix<CT>& , const GenVector<CT>& ,
        MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

#undef CT
}

#endif
