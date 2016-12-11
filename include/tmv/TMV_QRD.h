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
// This file contains the code for doing division using
// QR Decomposition.
//
// The basic idea of an QR decomposition is that any
// matrix A can be decomposed into a unitary matrix
// times an upper triangle matrix.
//
// A = Q R
//
// We do this using Householder transformations, which can
// be stored in a lower triangle matrix, thus other than the
// diagonal, which they both need, we can store them in the
// place of the original matrix. It is more convenient to
// keep the diagonal of Q in place and take the diagonal of R
// separately.
//
// If R is not singular, the solution to A x = b is found from
// Q R x = b
// R x = Qt b
// which can be solved by back-substitutaion.
//
// If m > n, this does not actually give a solution to A x = b,
// since Q Qt != I  (Q is only column orthogonal if m > n.)
// But it does give the value of x which minimizes the 2-norm
// of (A x - b)
//
// The 2-norm of v is the square root of vt v, so
// |A x - b|^2 =
// (A x - b)t (A x - b) =
// (xt At - bt) (A x - b) =
// xt At A x - bt A x - xt At b + bt b =
// xt Rt Qt Q R x - bt Q R x - xt Rt Qt b + bt b =
// |R x - Qt b|^2 + |b|^2 - |Qt b|^2
// Clearly the x which minimizes this is the solution of R x = Qt b.
//
// If R is singular, then you need QRP Decomposition (see TMV_QRPDiv.h).
//


#ifndef TMV_QRD_H
#define TMV_QRD_H

#include "tmv/TMV_Divider.h"
#include "tmv/TMV_BaseTriMatrix.h"

namespace tmv {

    // Decompose A (input as Q) into Q R
    // where Q is unitary and R is upper triangular.
    template <typename T>
    void QR_Decompose(MatrixView<T> Q, UpperTriMatrixView<T> R);

    // Decompose A into Q R, but don't return Q.
    // On output, R is returned as A.upperTri()
    template <typename T>
    void QR_Decompose(MatrixView<T> A);

    // Given that A0 = Q0 R0
    // Find R1, so that [ A0 ] = Q1 R1
    //                  [ A  ]
    // On input R is R0; on output it is R1.
    template <typename T>
    void QR_Update(UpperTriMatrixView<T> R, MatrixView<T> A);

    // This just adds a single row
    template <typename T>
    inline void QR_Update(UpperTriMatrixView<T> R, VectorView<T> z)
    { QR_Update(R,RowVectorViewOf(z)); }

    // The opposite of an Update:
    // Given that [ A0 ] = Q1 R1
    //            [ A  ]
    // Find R0, so that A0 = Q0 R0
    // On input R is R1; on output it is R0.
    template <typename T>
    void QR_Downdate(UpperTriMatrixView<T> R, MatrixView<T> A);

    // Remove a single row
    template <typename T>
    inline void QR_Downdate(UpperTriMatrixView<T> R, VectorView<T> z)
    { QR_Downdate(R,RowVectorViewOf(z)); }


    template <typename T, int A2>
    inline void QR_Decompose(MatrixView<T> Q, UpperTriMatrix<T,A2>& R)
    { QR_Decompose(Q,R.view()); }

    template <typename T, int A1>
    inline void QR_Decompose(Matrix<T,A1>& Q, UpperTriMatrixView<T> R)
    { QR_Decompose(Q.view(),R); }

    template <typename T, int A1, int A2>
    inline void QR_Decompose(Matrix<T,A1>& Q, UpperTriMatrix<T,A2>& R)
    { QR_Decompose(Q.view(),R.view()); }

    template <typename T, int A1>
    inline void QR_Decompose(Matrix<T,A1>& A)
    { QR_Decompose(A.view()); }

    template <typename T, int A2>
    inline void QR_Update(UpperTriMatrixView<T> R, Matrix<T,A2>& A)
    { QR_Update(R,A.view()); }

    template <typename T, int A1>
    inline void QR_Update(UpperTriMatrix<T,A1>& R, MatrixView<T> A)
    { QR_Update(R.view(),A); }

    template <typename T, int A1, int A2>
    inline void QR_Update(UpperTriMatrix<T,A1>& R, Matrix<T,A2>& A)
    { QR_Update(R.view(),A.view()); }

    template <typename T, int A2>
    inline void QR_Update(UpperTriMatrixView<T> R, Vector<T,A2>& A)
    { QR_Update(R,A.view()); }

    template <typename T, int A1>
    inline void QR_Update(UpperTriMatrix<T,A1>& R, VectorView<T> A)
    { QR_Update(R.view(),A); }

    template <typename T, int A1, int A2>
    inline void QR_Update(UpperTriMatrix<T,A1>& R, Vector<T,A2>& A)
    { QR_Update(R.view(),A.view()); }

    template <typename T, int A2>
    inline void QR_Downdate(UpperTriMatrixView<T> R, Matrix<T,A2>& A)
    { QR_Downdate(R,A.view()); }

    template <typename T, int A1>
    inline void QR_Downdate(UpperTriMatrix<T,A1>& R, MatrixView<T> A)
    { QR_Downdate(R.view(),A); }

    template <typename T, int A1, int A2>
    inline void QR_Downdate(UpperTriMatrix<T,A1>& R, Matrix<T,A2>& A)
    { QR_Downdate(R.view(),A.view()); }

    template <typename T, int A2>
    inline void QR_Downdate(UpperTriMatrixView<T> R, Vector<T,A2>& A)
    { QR_Downdate(R,A.view()); }

    template <typename T, int A1>
    inline void QR_Downdate(UpperTriMatrix<T,A1>& R, VectorView<T> A)
    { QR_Downdate(R.view(),A); }

    template <typename T, int A1, int A2>
    inline void QR_Downdate(UpperTriMatrix<T,A1>& R, Vector<T,A2>& A)
    { QR_Downdate(R.view(),A.view()); }

    template <typename T>
    class PackedQ;

    template <typename T>
    class QRDiv : public Divider<T>
    {

    public :

        //
        // Constructors
        //

        QRDiv(const GenMatrix<T>& A, bool _inplace);
        ~QRDiv();

        //
        // Div, DivEq
        //

        template <typename T1>
        void doLDivEq(MatrixView<T1> m) const;

        template <typename T1>
        void doRDivEq(MatrixView<T1> m) const;

        template <typename T1, typename T2>
        void doLDiv(const GenMatrix<T1>& m, MatrixView<T2> x) const;

        template <typename T1, typename T2>
        void doRDiv(const GenMatrix<T1>& m, MatrixView<T2> x) const;

#include "tmv/TMV_AuxAllDiv.h"

        //
        // Determinant, Inverse
        //

        T det() const;
        TMV_RealType(T) logDet(T* sign) const;
        template <typename T1>
        void doMakeInverse(MatrixView<T1> minv) const;
        void doMakeInverseATA(MatrixView<T> minv) const;
        bool isSingular() const;

        //
        // Access Decomposition
        //

        bool isTrans() const;
        PackedQ<T> getQ() const;
        ConstUpperTriMatrixView<T> getR() const;
        const GenMatrix<T>& getQR() const;
        const GenVector<T>& getBeta() const;

        bool checkDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

    protected :

        struct QRDiv_Impl;
        auto_ptr<QRDiv_Impl> pimpl;

        ptrdiff_t colsize() const;
        ptrdiff_t rowsize() const;

    private :

        QRDiv(const QRDiv<T>&);
        QRDiv<T>& operator=(const QRDiv<T>&);

    };

} // namespace mv


#endif
