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
// LU Decomposition.
//
// The basics of LU decomposition for band matrices are the same as
// for regular matrices.  However, there are a few wrinkles about doing
// it efficiently.
//
// We leave the details to the comments in TMV_BandLUDiv.cpp, but
// the main difference for the routines in this file is that L can
// be stored in a lower band matrix with m.nlo() subdiagonals.
// However, U needs m.nlo() + m.nhi() superdiagonals for its storage.
//
//


#ifndef TMV_BandLUD_H
#define TMV_BandLUD_H

#include "tmv/TMV_Divider.h"
#include "tmv/TMV_BaseBandMatrix.h"
#include "tmv/TMV_BaseMatrix.h"
#include "tmv/TMV_BaseTriMatrix.h"

namespace tmv {

    // Decompose A into P L U
    // L, U, and P must have the same size as A.
    // L should be UnitDiag.
    // U should be NonUnitDiag.
    template <typename T>
    void LU_Decompose(
        const GenBandMatrix<T>& A, LowerTriMatrixView<T> L,
        BandMatrixView<T> U, ptrdiff_t* P);

    // Do the decomposition in compressed form.
    template <typename T>
    void LU_Decompose(BandMatrixView<T> LUx, ptrdiff_t* P, ptrdiff_t Anhi);

    class Permutation;

    template <typename T>
    void LU_Decompose(
        const GenBandMatrix<T>& A, LowerTriMatrixView<T> L,
        BandMatrixView<T> U, Permutation& P);

    template <typename T>
    void LU_Decompose(BandMatrixView<T> LUx, Permutation& P, ptrdiff_t Anhi);

    template <typename T, int A1>
    inline void LU_Decompose(
        const GenBandMatrix<T>& A, LowerTriMatrix<T,A1>& L,
        BandMatrixView<T> U, Permutation& P)
    { LU_Decompose(A,L.view(),U,P); }

    template <typename T, int A2>
    inline void LU_Decompose(
        const GenBandMatrix<T>& A, LowerTriMatrixView<T> L,
        BandMatrix<T,A2>& U, Permutation& P)
    { LU_Decompose(A,L,U.view(),P); }

    template <typename T, int A1, int A2>
    inline void LU_Decompose(
        const GenBandMatrix<T>& A, LowerTriMatrix<T,A1>& L,
        BandMatrix<T,A2>& U, Permutation& P)
    { LU_Decompose(A,L.view(),U.view(),P); }

    template <typename T, int A1>
    inline void LU_Decompose(BandMatrix<T,A1>& LUx, Permutation& P, ptrdiff_t Anhi)
    { LU_Decompose(LUx.view(),P,Anhi); }

    template <typename T>
    class BandLUDiv : public Divider<T>
    {

    public :

        //
        // Constructors
        //

        BandLUDiv(const GenBandMatrix<T>& A, bool _inplace);
        BandLUDiv(const AssignableToBandMatrix<T>& A);
        ~BandLUDiv();

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
        LowerTriMatrix<T,UnitDiag> getL() const;
        ConstBandMatrixView<T> getU() const;
        const GenBandMatrix<T>& getLU() const;
        const Permutation& getP() const;

        bool checkDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

    private :

        struct BandLUDiv_Impl;
        auto_ptr<BandLUDiv_Impl> pimpl;

        ptrdiff_t colsize() const;
        ptrdiff_t rowsize() const;

    private :

        BandLUDiv(const BandLUDiv<T>&);
        BandLUDiv<T>& operator=(const BandLUDiv<T>&);

    };

} // namespace mv

#endif
