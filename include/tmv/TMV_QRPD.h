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
// QRP Decomposition.
//
// In a QR Decomposition, if R is singular, then back-substitution will
// fail in solving R x = Qt b.  However, there is a trick that still
// gives us a valid least squares solution to A x = b.
// This is to use column pivoting to obtain a decomposition:
// A = Q [ R11 R12 ] P
//       [  0   0  ]
// where R11 is upper triangular.
//
// With this decomposition,
// Q R P x = b
// R P x = Qt b
// Let z = P x (ie. we will solve for z first, then x = z / P = P * z)
// and c = Qt b
// Then R z = c
// Since R is singular, there is no exact solution, but for the least-squares
// problem, we just want to minimize |R z - c|.
// If z = [ z1 ] and c = [ c1 ] (with the same dimensional split as R11 and R12)
//        [ z2 ]         [ c2 ]
// then R z - c = [ R11 z1 + R12 z2 - c1 ]
//                [        -c2           ]
// The minimum norm will occur for any z2 if z1 = R11^-1 (c1 - R12 z2)
// So a solution can be found with z2 = 0.
// This solution then has minimum |A x - b| (ie. the least squares
// solution).
//
// For dividing from the other side with m > n, we can always find a
// possible solution (assuming R is non-singular):
//
// xt Q R P = bt
// xt Q R = bt Pt
// xt Q = bt Pt R^-1  (This is done by front-substitution)
// xt = bt Pt R^-1 Qt
//
// This solution does satisfy the equation xt A = bt, but
// it will not be the solution with minimum |x|.
// Use SVD for the min |x| solution.
//
// If R is singular, we of course have a problem with the front-substitution
// step.  This time, the QRP decomposition does not work quite as well.
// Take the front-substitution step (the Q and P steps don't affect
// this minimization), and write R as above:
//
// zt R = ct
// [ z1t z2t ] [ R11  R12 ] = [ c1t c2t ]
//             [  0    0  ]
// [ (z1t R11) (z1t R12) ] = [ c1t c2t]
// |zt R - ct|^2 = |z1t R11 - c1t|^2 + |z1t R12 - c2t|^2
// We can set z2t = 0, since it is arbitrary.
// But it is not so clear what the correct z1t is in this case.
// We take z1t = c1t R11^-1, as an approximate solution, but point out
// that this is not really correct. You should use SVD instead for
// the correct least squares solution in this case.
//
// You can control whether QRP does a strict reordering so that the
// diagonal elements of R are in decreasing order of absolute value
// (which I call Strict QRP), or whether they are just reordered well
// enough to put the correct zeros at the end (which I call Loose QRP) using
// the function tmv::UseStrictQRP();
// To turn the strict algorithm back off, use tmv::UseStrictQRP(false);
// The Loose QRP algorithm is faster and seems to be good enough for
// all the singular matrices I've tried it on, but someone might have a
// need for the strict behavior.
//


#ifndef TMV_QRPD_H
#define TMV_QRPD_H

#include "tmv/TMV_Divider.h"
#include "tmv/TMV_BaseTriMatrix.h"
#include "tmv/TMV_QRD.h"
#include "tmv/TMV_Permutation.h"

namespace tmv {

    // Decompose A (input as Q) into Q R P.
    template <typename T>
    void QRP_Decompose(
        MatrixView<T> Q, UpperTriMatrixView<T> R, ptrdiff_t* P, bool strict);

    template <typename T>
    void QRP_Decompose(
        MatrixView<T> QRx, VectorView<T> beta, ptrdiff_t* P, T& signdet, bool strict);

    class Permutation;

    // Default value of strict=false is given in TMV_Permutation.h
    // for these next two functions.
    template <typename T>
    void QRP_Decompose(
        MatrixView<T> Q, UpperTriMatrixView<T> R, Permutation& P, bool strict);

    template <typename T>
    void QRP_Decompose(
        MatrixView<T> QRx, VectorView<T> beta, Permutation& P, T& signdet,
        bool strict);

    // Decompose A into Q R P, but don't return Q or P.
    // R is returned as A.upperTri().
    template <typename T>
    void QRP_Decompose(MatrixView<T> A, bool strict=false);

    template <typename T, int A2>
    inline void QRP_Decompose(
        MatrixView<T> Q, UpperTriMatrix<T,A2>& R, Permutation& P,
        bool strict=false)
    { QRP_Decompose(Q,R.view(),P,strict); }

    template <typename T, int A1>
    inline void QRP_Decompose(
        Matrix<T,A1>& Q, UpperTriMatrixView<T> R, Permutation& P,
        bool strict=false)
    { QRP_Decompose(Q.view(),R,P,strict); }

    template <typename T, int A1, int A2>
    inline void QRP_Decompose(
        Matrix<T,A1>& Q, UpperTriMatrix<T,A2>& R, Permutation& P,
        bool strict=false)
    { QRP_Decompose(Q.view(),R.view(),P,strict); }

    template <typename T, int A1>
    inline void QRP_Decompose(Matrix<T,A1>& A, bool strict=false)
    { QRP_Decompose(A.view(),strict); }

    struct QRP_StrictSingleton
    {
        // Technically, I think this isn't thread safe, but I'd be
        // pretty shocked if people were having multiple threads
        // call this funtion at the same time.
        static inline bool& inst()
        {
            static bool strict;
            return strict;
        }
    };

    inline void UseStrictQRP(bool newstrict=true)
    { QRP_StrictSingleton::inst() = newstrict; }

    inline bool QRP_IsStrict()
    { return QRP_StrictSingleton::inst(); }

    template <typename T>
    class QRPDiv : public Divider<T>
    {

    public :

        //
        // Constructors
        //

        QRPDiv(const GenMatrix<T>& A, bool _inplace);
        ~QRPDiv();

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
        PackedQ<T> getQ() const;
        ConstUpperTriMatrixView<T> getR() const;
        const GenMatrix<T>& getQRx() const;
        const GenVector<T>& getBeta() const;
        const Permutation& getP() const;

        bool checkDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

    protected :

        struct QRPDiv_Impl;
        auto_ptr<QRPDiv_Impl> pimpl;

        ptrdiff_t colsize() const;
        ptrdiff_t rowsize() const;

    private :

        QRPDiv(const QRPDiv<T>&);
        QRPDiv<T>& operator=(const QRPDiv<T>&);

    };

} // namespace mv

#endif
