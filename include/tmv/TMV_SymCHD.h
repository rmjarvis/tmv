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


#ifndef TMV_SymCHD_H
#define TMV_SymCHD_H

#include "tmv/TMV_SymDivider.h"
#include "tmv/TMV_BaseTriMatrix.h"
#include "tmv/TMV_BaseBandMatrix.h"

namespace tmv {

    // Decompose A into L*Lt
    // On output L = A.lowerTri()
    template <class T> 
    void CH_Decompose(SymMatrixView<T> A);

    // Decompose A into U*P
    // where U is unitary and P is positive definite
    // On ouput U = A
    template <class T> 
    void PolarDecompose(MatrixView<T> A, SymMatrixView<T> P);

    template <class T> 
    void PolarDecompose(
        const GenBandMatrix<T>& A, MatrixView<T> U, SymMatrixView<T> P);

    template <class T> 
    void SquareRoot(SymMatrixView<T> A);

    template <class T, int A1> 
    inline void CH_Decompose(HermMatrix<T,A1>& A)
    { CH_Decompose(A.view()); }

    template <class T, int A1> 
    inline void CH_Decompose(SymMatrix<T,A1>& A)
    { CH_Decompose(A.view()); }

    template <class T, int A2> 
    inline void PolarDecompose(MatrixView<T> A, HermMatrix<T,A2>& P)
    { PolarDecompose(A,P.view()); }

    template <class T, int A2> 
    inline void PolarDecompose(MatrixView<T> A, SymMatrix<T,A2>& P)
    { PolarDecompose(A,P.view()); }

    template <class T, int A1> 
    inline void PolarDecompose(Matrix<T,A1>& A, SymMatrixView<T> P)
    { PolarDecompose(A.view(),P); }

    template <class T, int A1, int A2> 
    inline void PolarDecompose(Matrix<T,A1>& A, HermMatrix<T,A2>& P)
    { PolarDecompose(A.view(),P.view()); }

    template <class T, int A1, int A2> 
    inline void PolarDecompose(Matrix<T,A1>& A, SymMatrix<T,A2>& P)
    { PolarDecompose(A.view(),P.view()); }

    template <class T, int A2> 
    inline void PolarDecompose(
        const GenBandMatrix<T>& A, MatrixView<T> U, HermMatrix<T,A2>& P)
    { PolarDecompose(A,U,P.view()); }

    template <class T, int A2> 
    inline void PolarDecompose(
        const GenBandMatrix<T>& A, MatrixView<T> U, SymMatrix<T,A2>& P)
    { PolarDecompose(A,U,P.view()); }

    template <class T, int A1> 
    inline void PolarDecompose(
        const GenBandMatrix<T>& A, Matrix<T,A1>& U, SymMatrixView<T> P)
    { PolarDecompose(A,U.view(),P); }

    template <class T, int A1, int A2> 
    inline void PolarDecompose(
        const GenBandMatrix<T>& A, Matrix<T,A1>& U, HermMatrix<T,A2>& P)
    { PolarDecompose(A,U.view(),P.view()); }

    template <class T, int A1, int A2> 
    inline void PolarDecompose(
        const GenBandMatrix<T>& A, Matrix<T,A1>& U, SymMatrix<T,A2>& P)
    { PolarDecompose(A,U.view(),P.view()); }

    template <class T, int A1> 
    inline void SquareRoot(SymMatrix<T,A1>& A)
    { SquareRoot(A.view()); }

    template <class T, int A1> 
    inline void SquareRoot(HermMatrix<T,A1>& A)
    { SquareRoot(A.view()); }


    template <class T> 
    class HermCHDiv : public SymDivider<T> 
    {

    public :

        //
        // Constructors
        //

        HermCHDiv(const GenSymMatrix<T>& A, bool _inplace);
        ~HermCHDiv();

        //
        // Divider Versions of DivEq and Div
        //

        template <class T1> 
        void doLDivEq(MatrixView<T1> m) const;

        template <class T1> 
        void doRDivEq(MatrixView<T1> m) const;

        template <class T1, class T2> 
        void doLDiv(const GenMatrix<T1>& m1, MatrixView<T2> m0) const;

        template <class T1, class T2> 
        void doRDiv(const GenMatrix<T1>& m1, MatrixView<T2> m0) const;

        //
        // Determinant, Inverse
        //

        T det() const;
        TMV_RealType(T) logDet(T* sign) const;
        template <class T1> 
        void doMakeInverse(MatrixView<T1> minv) const;

        template <class T1> 
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

        const ConstLowerTriMatrixView<T> getL() const;
        const GenSymMatrix<T>& getLL() const;

        bool checkDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

    private :
        struct HermCHDiv_Impl;
        std::auto_ptr<HermCHDiv_Impl> pimpl;

        ptrdiff_t colsize() const;
        ptrdiff_t rowsize() const;

    private :

        HermCHDiv(const HermCHDiv<T>&);
        HermCHDiv<T>& operator=(const HermCHDiv<T>&);

    }; // HermCHDiv

} // namespace tmv

#endif
