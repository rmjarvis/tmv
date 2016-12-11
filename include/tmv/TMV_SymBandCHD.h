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
// Cholesky Decomposition.
//
// The algorithm is much like the LU decomposition, but we don't do
// any pivoting, and since the source matrix is symmetric, L = LT
// (or for Hermition, L = Lt).
//


#ifndef TMV_SymBandCHD_H
#define TMV_SymBandCHD_H

#include "tmv/TMV_SymDivider.h"
#include "tmv/TMV_BaseSymBandMatrix.h"
#include "tmv/TMV_BaseDiagMatrix.h"

namespace tmv {

    // Decompose A into L*Lt.
    // On output, L = A.lowerBand().
    template <typename T>
    void CH_Decompose(SymBandMatrixView<T> A);

    // Decompose Tridiagonal A into L D Lt
    template <typename T>
    void LDL_Decompose(SymBandMatrixView<T> A);

    template <typename T, int A1>
    inline void CH_Decompose(HermBandMatrix<T,A1>& A)
    { CH_Decompose(A.view()); }

    template <typename T, int A1>
    inline void CH_Decompose(SymBandMatrix<T,A1>& A)
    { CH_Decompose(A.view()); }

    template <typename T, int A1>
    inline void LDL_Decompose(HermBandMatrix<T,A1>& A)
    { LDL_Decompose(A.view()); }

    template <typename T, int A1>
    inline void LDL_Decompose(SymBandMatrix<T,A1>& A)
    { LDL_Decompose(A.view()); }


    template <typename T>
    class HermBandCHDiv : public SymDivider<T>
    {

    public :

        //
        // Constructors
        //

        HermBandCHDiv(const GenSymBandMatrix<T>& A, bool _inplace);
        ~HermBandCHDiv();

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

        //
        // Determinant, Inverse
        //

        T det() const;
        TMV_RealType(T) logDet(T* sign) const;
        template <typename T1>
        void doMakeInverse(MatrixView<T1> minv) const;
        template <typename T1>
        void doMakeInverse(SymMatrixView<T1> minv) const;
        inline void makeInverse(SymMatrixView<TMV_RealType(T)> sinv) const
        {
            TMVAssert(isReal(T()));
            TMVAssert(sinv.size() == colsize());
            doMakeInverse(sinv);
        }
        inline void makeInverse(SymMatrixView<TMV_ComplexType(T)> sinv) const
        {
            TMVAssert(sinv.size() == colsize());
            TMVAssert(sinv.isherm());
            doMakeInverse(sinv);
        }
        void doMakeInverseATA(MatrixView<T> minv) const;
        bool isSingular() const;

#include "tmv/TMV_AuxAllDiv.h"

        //
        // Access Decomposition
        //

        const BandMatrix<T> getL() const;
        const DiagMatrix<T> getD() const;
        const GenSymBandMatrix<T>& getLL() const;

        bool checkDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

    private :
        struct HermBandCHDiv_Impl;
        auto_ptr<HermBandCHDiv_Impl> pimpl;

        ptrdiff_t colsize() const;
        ptrdiff_t rowsize() const;

    private :

        HermBandCHDiv(const HermBandCHDiv<T>&);
        HermBandCHDiv<T>& operator=(const HermBandCHDiv<T>&);

    };

} // namespace tmv

#endif
