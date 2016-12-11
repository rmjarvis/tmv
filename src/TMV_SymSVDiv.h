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

#ifndef TMV_SymSVDiv_H
#define TMV_SymSVDiv_H

#include "tmv/TMV_BaseSymMatrix.h"
#include "TMV_SVDiv.h"

namespace tmv {

#define RT TMV_RealType(T)
#define RT1 TMV_RealType(T1)

    template <typename T, typename Td> 
    void Tridiagonalize(
        SymMatrixView<T> A, VectorView<T> beta,
        VectorView<Td> D, VectorView<RT> E, T& signdet);

    template <typename T> 
    void EigenFromTridiagonal(
        MatrixView<T> U, VectorView<RT> D, VectorView<RT> E);

    template <typename T> 
    void UnsortedHermEigen(MatrixView<T> U, VectorView<RT> S);
    template <typename T> 
    void UnsortedEigen(SymMatrixView<T> A, VectorView<RT> S);

    template <typename T> 
    void HermSV_Decompose(MatrixView<T> U, DiagMatrixView<RT> S);
    template <typename T> 
    void SymSV_Decompose(
        MatrixView<T> U,
        DiagMatrixView<RT> S, MatrixView<T> Vt, RT& logdet, T& signdet);
    template <typename T> 
    void SV_Decompose(
        SymMatrixView<T> A, DiagMatrixView<RT> S, MatrixView<T> Vt);

    template <typename T, typename T1> 
    void HermSV_Inverse(
        const GenMatrix<T1>& U, const GenDiagMatrix<RT1>& S, 
        ptrdiff_t kmax, SymMatrixView<T> sinv);
    template <typename T, typename T1> 
    void SymSV_Inverse(
        const GenMatrix<T1>& U, const GenDiagMatrix<RT1>& S,
        const GenMatrix<T1>& Vt, ptrdiff_t kmax, SymMatrixView<T> sinv);

    template <typename T> 
    void HermTridiagonalChopSmallElements(
        VectorView<T> D, VectorView<T> E);

    template <typename T> 
    void EigenFromTridiagonal_QR(
        MatrixView<T> U, VectorView<RT> D, VectorView<RT> E);

    template <typename T> 
    void EigenFromTridiagonal_DC(
        MatrixView<T> U, VectorView<RT> D, VectorView<RT> E, bool UisI);

    template <typename T> 
    void FindDCEigenValues(
        Vector<T>& S, const T rho, const GenVector<T>& D, const GenVector<T>& z,
        Matrix<T,ColMajor>& diffmat);

    template <typename T> 
    void FindDCEigenValues(
        Vector<T>& S, const T rho, const GenVector<T>& D, const GenVector<T>& z);

#undef RT

#define CT std::complex<T>

    // Again, Microsoft Visual C++ needs this extra parameter to 
    // get the resolution of the overload right.
    template <typename T, typename T1> 
    inline void CallHermSV_Inverse(
        T , const GenMatrix<T1>& U, const GenDiagMatrix<RT1>& S,
        ptrdiff_t kmax, SymMatrixView<T> sinv)
    { HermSV_Inverse(U,S,kmax,sinv); }

    template <typename T, typename T1> 
    inline void CallHermSV_Inverse(
        T , const GenMatrix<T1>& U, const GenDiagMatrix<RT1>& S, 
        ptrdiff_t kmax, SymMatrixView<CT> sinv)
    { HermSV_Inverse(U,S,kmax,sinv); }

    template <typename T, typename T1> 
    inline void CallSymSV_Inverse(
        T , const GenMatrix<T1>& U, const GenDiagMatrix<RT1>& S,
        const GenMatrix<T1>& Vt, ptrdiff_t kmax, SymMatrixView<T> sinv)
    { SymSV_Inverse(U,S,Vt,kmax,sinv); }

    template <typename T, typename T1> 
    inline void CallSymSV_Inverse(
        T , const GenMatrix<T1>& U, const GenDiagMatrix<RT1>& S,
        const GenMatrix<T1>& Vt, ptrdiff_t kmax, SymMatrixView<CT> sinv)
    { SymSV_Inverse(U,S,Vt,kmax,sinv); }

    // Specialize disallowed complex combinations:
    template <typename T>
    inline void CallHermSV_Inverse(
        CT , const GenMatrix<CT>& , const GenDiagMatrix<T>& , 
        ptrdiff_t , SymMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T>
    inline void CallSymSV_Inverse(
        CT , const GenMatrix<CT>& , const GenDiagMatrix<T>& ,
        const GenMatrix<CT>& , ptrdiff_t , SymMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

#undef CT
 
}

#endif
