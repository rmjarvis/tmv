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


//---------------------------------------------------------------------------
//
// This file contains the code for doing division of BandMatrices using
// QR Decomposition.
//
// I don't implement the QRP method, since the permutations screw up the
// band structure.  It could be implemented using a similar technique
// as I used for the BandLUDiv class, but for now if you want to use it
// you need to copy the band matrix to a regular matrix first.
//


#ifndef TMV_BandQRD_H
#define TMV_BandQRD_H

#include "tmv/TMV_Divider.h"
#include "tmv/TMV_BaseBandMatrix.h"
#include "tmv/TMV_BaseMatrix.h"

namespace tmv {

    // Decompose A into Q R
    // where Q is unitary and R is upper banded.
    // R must have R.nhi() >= A.nlo()+A.nhi() and R.nlo() >= 0
    // Also, A must have A.nrows() >= A.ncols()
    template <typename T>
    void QR_Decompose(
        const GenBandMatrix<T>& A, MatrixView<T> Q, BandMatrixView<T> R);

    // The same, but don't return Q
    template <typename T>
    void QR_Decompose(const GenBandMatrix<T>& A, BandMatrixView<T> R);

    template <typename T, int A2>
    inline void QR_Decompose(
        const GenBandMatrix<T>& A, MatrixView<T> Q, BandMatrix<T,A2>& R)
    { QR_Decompose(A,Q,R.view()); }

    template <typename T, int A1>
    inline void QR_Decompose(
        const GenBandMatrix<T>& A, Matrix<T,A1>& Q, BandMatrixView<T> R)
    { QR_Decompose(A,Q.view(),R); }

    template <typename T, int A1, int A2>
    inline void QR_Decompose(
        const GenBandMatrix<T>& A, Matrix<T,A1>& Q, BandMatrix<T,A2>& R)
    { QR_Decompose(A,Q.view(),R.view()); }

    template <typename T, int A2>
    inline void QR_Decompose(const GenBandMatrix<T>& A, BandMatrix<T,A2>& R)
    { QR_Decompose(A,R.view()); }

    template <typename T>
    class BandQRDiv : public Divider<T>
    {

    public :

        //
        // Constructors
        //

        BandQRDiv(const GenBandMatrix<T>& A, bool _inplace);
        ~BandQRDiv();

        //
        // Div, DivEq
        //

        template <typename T1>
        void doLDivEq(MatrixView<T1> m) const;
        template <typename T1>
        void doRDivEq(MatrixView<T1> m) const;
        template <typename T1, typename T2>
        void doLDiv(
            const GenMatrix<T1>& m, MatrixView<T2> x) const;
        template <typename T1, typename T2>
        void doRDiv(
            const GenMatrix<T1>& m, MatrixView<T2> x) const;

        //
        // Determinant, Inverse
        //

        T det() const;
        TMV_RealType(T) logDet(T* sign) const;
        template <typename T1>
        void doMakeInverse(MatrixView<T1> minv) const;
        void doMakeInverseATA(MatrixView<T> minv) const;
        bool isSingular() const;

#include "tmv/TMV_AuxAllDiv.h"

        //
        // Access Decomposition
        //

        bool isTrans() const;
        Matrix<T> getQ() const;
        ConstBandMatrixView<T> getR() const;
        const GenBandMatrix<T>& getQR() const;
        const GenVector<T>& getQBeta() const;

        bool checkDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

    private :

        struct BandQRDiv_Impl;
        auto_ptr<BandQRDiv_Impl> pimpl;

        ptrdiff_t colsize() const;
        ptrdiff_t rowsize() const;

    private :

        BandQRDiv(const BandQRDiv<T>&);
        BandQRDiv<T>& operator=(const BandQRDiv<T>&);

    };

} // namespace mv

#endif
