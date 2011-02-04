///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2009                                                 //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file LICENSE.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#ifndef TMV_BandQRDiv_H
#define TMV_BandQRDiv_H

#include "tmv/TMV_BaseBandMatrix.h"

namespace tmv {

    template <class T> 
    void QR_Decompose(
        const BandMatrixView<T>& QRx, const VectorView<T>& Qbeta, T& signdet);

    template <class T> 
    void GetQFromBandQR(
        const MatrixView<T>& QRx, const GenVector<T>& Qbeta, int nlo);

    template <class T, class T1> 
    void QR_LDivEq(
        const GenBandMatrix<T1>& QRx, const GenVector<T1>& Qbeta,
        const MatrixView<T>& m);
    template <class T, class T1> 
    void QR_RDivEq(
        const GenBandMatrix<T1>& QRx, const GenVector<T1>& Qbeta,
        const MatrixView<T>& m);

    template <class T, class T1, class T2> 
    void QR_LDiv(
        const GenBandMatrix<T1>& QR, const GenVector<T1>& Qbeta,
        const GenMatrix<T2>& m, const MatrixView<T>& x);
    template <class T, class T1, class T2> 
    void QR_RDiv(
        const GenBandMatrix<T1>& QR, const GenVector<T1>& Qbeta,
        const GenMatrix<T2>& m, const MatrixView<T>& x);

    template <class T, class T1> 
    void Q_LDivEq(
        const GenBandMatrix<T1>& Q, const GenVector<T1>& Qbeta,
        const MatrixView<T>& m);
    template <class T, class T1> 
    void Q_RDivEq(
        const GenBandMatrix<T1>& Q, const GenVector<T1>& Qbeta,
        const MatrixView<T>& m);

    template <class T, class T1> 
    void QR_Inverse(
        const GenBandMatrix<T1>& QRx, const GenVector<T1>& Qbeta,
        const MatrixView<T>& m);


    // Specialize disallowed complex combinations:
#define CT std::complex<T>

    template <class T>
    inline void QR_LDivEq(
        const GenBandMatrix<CT>& QRx, const GenVector<CT>& Qbeta,
        const MatrixView<T>& m)
    { TMVAssert(TMV_FALSE); }
    template <class T>
    inline void QR_RDivEq(
        const GenBandMatrix<CT>& QRx, const GenVector<CT>& Qbeta,
        const MatrixView<T>& m)
    { TMVAssert(TMV_FALSE); }

    template <class T>
    inline void QR_LDiv(
        const GenBandMatrix<CT>& QR, const GenVector<CT>& Qbeta,
        const GenMatrix<CT>& m, const MatrixView<T>& x)
    { TMVAssert(TMV_FALSE); }
    template <class T>
    inline void QR_LDiv(
        const GenBandMatrix<CT>& QR, const GenVector<CT>& Qbeta,
        const GenMatrix<T>& m, const MatrixView<T>& x)
    { TMVAssert(TMV_FALSE); }
    template <class T>
    inline void QR_LDiv(
        const GenBandMatrix<T>& QR, const GenVector<T>& Qbeta,
        const GenMatrix<CT>& m, const MatrixView<T>& x)
    { TMVAssert(TMV_FALSE); }
    template <class T>
    inline void QR_RDiv(
        const GenBandMatrix<CT>& QR, const GenVector<CT>& Qbeta,
        const GenMatrix<CT>& m, const MatrixView<T>& x)
    { TMVAssert(TMV_FALSE); }
    template <class T>
    inline void QR_RDiv(
        const GenBandMatrix<CT>& QR, const GenVector<CT>& Qbeta,
        const GenMatrix<T>& m, const MatrixView<T>& x)
    { TMVAssert(TMV_FALSE); }
    template <class T>
    inline void QR_RDiv(
        const GenBandMatrix<T>& QR, const GenVector<T>& Qbeta,
        const GenMatrix<CT>& m, const MatrixView<T>& x)
    { TMVAssert(TMV_FALSE); }

    template <class T>
    inline void Q_LDivEq(
        const GenBandMatrix<CT>& Q, const GenVector<CT>& Qbeta,
        const MatrixView<T>& m)
    { TMVAssert(TMV_FALSE); }
    template <class T>
    inline void Q_RDivEq(
        const GenBandMatrix<CT>& Q, const GenVector<CT>& Qbeta,
        const MatrixView<T>& m)
    { TMVAssert(TMV_FALSE); }

    template <class T>
    inline void QR_Inverse(
        const GenBandMatrix<CT>& QRx, const GenVector<CT>& Qbeta,
        const MatrixView<T>& m)
    { TMVAssert(TMV_FALSE); }

#undef CT
}

#endif
