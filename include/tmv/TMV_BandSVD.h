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
#ifndef TMV_BandSVD_H
#define TMV_BandSVD_H

#include "tmv/TMV_Divider.h"
#include "tmv/TMV_BaseBandMatrix.h"
#include "tmv/TMV_BaseDiagMatrix.h"

namespace tmv {

    // Decompose A into U S Vt
    // where S is a diagonal real matrix, and U,Vt are unitary matrices.
    // U,S,Vt are N x N
    template <typename T>
    void SV_Decompose(
        const GenBandMatrix<T>& A,
        MatrixView<T> U, DiagMatrixView<TMV_RealType(T)> S, MatrixView<T> Vt);

    // The same decomposition, but don't return U and/or Vt
    template <typename T>
    void SV_Decompose(
        const GenBandMatrix<T>& A,
        MatrixView<T> U, DiagMatrixView<TMV_RealType(T)> S);
    template <typename T>
    void SV_Decompose(
        const GenBandMatrix<T>& A,
        DiagMatrixView<TMV_RealType(T)> S, MatrixView<T> Vt);
    template <typename T>
    void SV_Decompose(
        const GenBandMatrix<T>& A, DiagMatrixView<TMV_RealType(T)> S);

    template <typename T, int A3>
    inline void SV_Decompose(
        const GenBandMatrix<T>& A,
        MatrixView<T> U, DiagMatrixView<TMV_RealType(T)> S, Matrix<T,A3>& Vt)
    { SV_Decompose(A,U,S,Vt.view()); }

    template <typename T, int A2>
    inline void SV_Decompose(
        const GenBandMatrix<T>& A,
        MatrixView<T> U, DiagMatrix<TMV_RealType(T),A2>& S, MatrixView<T> Vt)
    { SV_Decompose(A,U,S.view(),Vt); }

    template <typename T, int A2, int A3>
    inline void SV_Decompose(
        const GenBandMatrix<T>& A,
        MatrixView<T> U, DiagMatrix<TMV_RealType(T),A2>& S, Matrix<T,A3>& Vt)
    { SV_Decompose(A,U,S.view(),Vt.view()); }

    template <typename T, int A1>
    inline void SV_Decompose(
        const GenBandMatrix<T>& A,
        Matrix<T,A1>& U, DiagMatrixView<TMV_RealType(T)> S, MatrixView<T> Vt)
    { SV_Decompose(A,U.view(),S,Vt); }

    template <typename T, int A1, int A3>
    inline void SV_Decompose(
        const GenBandMatrix<T>& A,
        Matrix<T,A1>& U, DiagMatrixView<TMV_RealType(T)> S, Matrix<T,A3>& Vt)
    { SV_Decompose(A,U.view(),S,Vt.view()); }

    template <typename T, int A1, int A2>
    inline void SV_Decompose(
        const GenBandMatrix<T>& A,
        Matrix<T,A1>& U, DiagMatrix<TMV_RealType(T),A2>& S, MatrixView<T> Vt)
    { SV_Decompose(A,U.view(),S.view(),Vt); }

    template <typename T, int A1, int A2, int A3>
    inline void SV_Decompose(
        const GenBandMatrix<T>& A,
        Matrix<T,A1>& U, DiagMatrix<TMV_RealType(T),A2>& S, Matrix<T,A3>& Vt)
    { SV_Decompose(A,U.view(),S.view(),Vt.view()); }

    template <typename T, int A2>
    inline void SV_Decompose(
        const GenBandMatrix<T>& A,
        MatrixView<T> U, DiagMatrix<TMV_RealType(T),A2>& S)
    { SV_Decompose(A,U,S.view()); }

    template <typename T, int A1>
    inline void SV_Decompose(
        const GenBandMatrix<T>& A,
        Matrix<T,A1>& U, DiagMatrixView<TMV_RealType(T)> S)
    { SV_Decompose(A,U.view(),S); }

    template <typename T, int A1, int A2>
    inline void SV_Decompose(
        const GenBandMatrix<T>& A,
        Matrix<T,A1>& U, DiagMatrix<TMV_RealType(T),A2>& S)
    { SV_Decompose(A,U.view(),S.view()); }

    template <typename T, int A3>
    inline void SV_Decompose(
        const GenBandMatrix<T>& A,
        DiagMatrixView<TMV_RealType(T)> S, Matrix<T,A3>& Vt)
    { SV_Decompose(A,S,Vt.view()); }

    template <typename T, int A2>
    inline void SV_Decompose(
        const GenBandMatrix<T>& A,
        DiagMatrix<TMV_RealType(T),A2>& S, MatrixView<T> Vt)
    { SV_Decompose(A,S.view(),Vt); }

    template <typename T, int A2, int A3>
    inline void SV_Decompose(
        const GenBandMatrix<T>& A,
        DiagMatrix<TMV_RealType(T),A2>& S, Matrix<T,A3>& Vt)
    { SV_Decompose(A,S.view(),Vt.view()); }

    template <typename T, int A2>
    inline void SV_Decompose(
        const GenBandMatrix<T>& A, DiagMatrix<TMV_RealType(T),A2>& S)
    { SV_Decompose(A,S.view()); }


    template <typename T>
    class BandSVDiv : public Divider<T>
    {

    public :

        //
        // Constructors
        //

        BandSVDiv(const GenBandMatrix<T>& A);
        ~BandSVDiv();

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
        TMV_RealType(T) norm2() const;
        TMV_RealType(T) condition() const;

        //
        // Determine which (if any) S values to zero out
        //

        void thresh(TMV_RealType(T) toler, std::ostream* debugout=0) const;
        void top(ptrdiff_t neigen, std::ostream* debugout=0) const;
        ptrdiff_t getKMax() const;

        //
        // Access Decomposition
        //

        ConstMatrixView<T> getU() const;
        ConstDiagMatrixView<TMV_RealType(T)> getS() const;
        ConstMatrixView<T> getVt() const;

        bool checkDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

    protected:

        struct BandSVDiv_Impl;
        auto_ptr<BandSVDiv_Impl> pimpl;

        ptrdiff_t colsize() const;
        ptrdiff_t rowsize() const;

    private :

        BandSVDiv(const BandSVDiv<T>&);
        BandSVDiv<T>& operator=(const BandSVDiv<T>&);

    }; // BandSVDiv

} // namespace tmv;

#endif
