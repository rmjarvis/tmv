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
#ifndef TMV_SymSVDiv_H
#define TMV_SymSVDiv_H

#include "tmv/TMV_BaseSymMatrix.h"
#include "TMV_SVDiv.h"

namespace tmv {

#define RT TMV_RealType(T)
#define RT1 TMV_RealType(T1)

    template <class T, class Td> 
    void Tridiagonalize(
        const SymMatrixView<T>& A, const VectorView<T>& beta,
        const VectorView<Td>& D, const VectorView<RT>& E, T& signdet);

    template <class T> 
    void EigenFromTridiagonal(
        MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E);

    template <class T> 
    void UnsortedHermEigen(const MatrixView<T>& U, const VectorView<RT>& S);
    template <class T> 
    void UnsortedEigen(const SymMatrixView<T>& A, const VectorView<RT>& S);

    template <class T> 
    void HermSV_Decompose(const MatrixView<T>& U, const DiagMatrixView<RT>& S);
    template <class T> 
    void SymSV_Decompose(
        const MatrixView<T>& U,
        const DiagMatrixView<RT>& S, MVP<T> V, RT& logdet, T& signdet);
    template <class T> 
    void SV_Decompose(
        const SymMatrixView<T>& A, const DiagMatrixView<RT>& S, MVP<T> V);

    template <class T, class T1> 
    void HermSV_Inverse(
        const GenMatrix<T1>& U, const GenDiagMatrix<RT1>& S, int kmax,
        const SymMatrixView<T>& sinv);
    template <class T, class T1> 
    void SymSV_Inverse(
        const GenMatrix<T1>& U, const GenDiagMatrix<RT1>& S,
        const GenMatrix<T1>& V, int kmax, const SymMatrixView<T>& sinv);

    template <class T> 
    void HermTridiagonalChopSmallElements(
        const VectorView<T>& D, const VectorView<T>& E);

    template <class T> 
    void EigenFromTridiagonal_QR(
        MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E);

    template <class T> 
    void EigenFromTridiagonal_DC(
        MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E, bool UisI);

    template <class T> 
    void FindDCEigenValues(
        Vector<T>& S, const T rho, const GenVector<T>& D, const GenVector<T>& z,
        Matrix<T,ColMajor>& diffmat);

    template <class T> 
    void FindDCEigenValues(
        Vector<T>& S, const T rho, const GenVector<T>& D, const GenVector<T>& z);

#undef RT

    // Specialize disallowed complex combinations:
#define CT std::complex<T>

    template <class T>
    void HermSV_Inverse(
        const GenMatrix<CT>& U, const GenDiagMatrix<T>& S, int kmax,
        const SymMatrixView<T>& sinv)
    { TMVAssert(TMV_FALSE); }
    template <class T>
    void SymSV_Inverse(
        const GenMatrix<CT>& U, const GenDiagMatrix<T>& S,
        const GenMatrix<CT>& V, int kmax, const SymMatrixView<T>& sinv)
    { TMVAssert(TMV_FALSE); }

#undef CT
 
}

#endif
