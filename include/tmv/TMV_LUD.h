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
// LU Decomposition.
//
// The name LU Decomposition is traditional, but somewhat
// annoying.  Usually U represents a unitary matrix, not an
// upper traiangular matrix.  The latter are usually represented
// with R.  (for "Right Triangular")
// For example, in a QR decomposition, the R is an upper
// trianular matrix.  (Q also typically represents unitary matrices.)
// However, I will use U rather than R here, since that is
// the usual representation in this context.
//
// The basic idea of an LU decomposition is that any
// square matrix A can be decomposed into a lower triangular
// matrix time an upper triangular matrix.
//
// For stability reasons, we actually decompose a permutation
// of A instead, so:
//
// A = P L U
//
// Only one of L or U needs a non-unit diagonal, so we choose L to
// have unit diagonal, and U to have the non-unit diagonal.
//
// This means that we can store L and U both in a square matrix
// the same size as A, with L being the elements below the diagonal
// and U being the elements above and including the diagonal.
//
// The determinant of A can be calculated easily from the LU
// decomposition:
//
// det(P) * det(A) = det(L) * det(U)
// +-1 * det(A) = 1 * det(U)
// As we calculate the decomposition, we keep track of whether
// det(P) is +-1
// The determinant of U is just the product of the diagonal elements.
// So the determinant of A is just det(P) times the diagonal elements
// of U.
//


#ifndef TMV_LUD_H
#define TMV_LUD_H

#include "tmv/TMV_Divider.h"
#include "tmv/TMV_BaseMatrix.h"
#include "tmv/TMV_BaseTriMatrix.h"

namespace tmv {

    // Decompose A into P * L * U
    // L is returned as A.lowerTri(UnitDiag).
    // U is returned as A.upperTri(NonUnitDiag).
    template <typename T>
    void LU_Decompose(MatrixView<T> A, ptrdiff_t* P);

    class Permutation;

    template <typename T>
    void LU_Decompose(MatrixView<T> A, Permutation& P);

    template <typename T, int A1>
    inline void LU_Decompose(Matrix<T,A1>& A, Permutation& P)
    { LU_Decompose(A.view(),P); }

    template <typename T>
    class LUDiv : public Divider<T>
    {

    public :

        //
        // Constructors
        //

        LUDiv(const GenMatrix<T>& A, bool _inplace);
        ~LUDiv();

        //
        // Divider Versions of DivEq and Div
        //

        template <typename T1>
        void doLDivEq(MatrixView<T1> m) const;

        template <typename T1>
        void doRDivEq(MatrixView<T1> m) const;

        template <typename T1, typename T2>
        void doLDiv(const GenMatrix<T1>& m1, MatrixView<T2> m0) const;

        template <typename T1, typename T2>
        void doRDiv(const GenMatrix<T1>& m1, MatrixView<T2> m0) const;

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
        ConstLowerTriMatrixView<T> getL() const;
        ConstUpperTriMatrixView<T> getU() const;
        const GenMatrix<T>& getLU() const;
        const Permutation& getP() const;

        bool checkDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

    private :

        struct LUDiv_Impl;
        auto_ptr<LUDiv_Impl> pimpl;

        ptrdiff_t colsize() const;
        ptrdiff_t rowsize() const;

    private :

        LUDiv(const LUDiv<T>&);
        LUDiv<T>& operator=(const LUDiv<T>&);

    };

} // namespace tmv

#endif
