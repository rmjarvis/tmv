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
        SymMatrixView<T> A, VectorView<T> beta,
        VectorView<Td> D, VectorView<RT> E, T& signdet);

    template <class T> 
    void EigenFromTridiagonal(
        MVP<T> U, VectorView<RT> D, VectorView<RT> E);

    template <class T> 
    void UnsortedHermEigen(MatrixView<T> U, VectorView<RT> S);
    template <class T> 
    void UnsortedEigen(SymMatrixView<T> A, VectorView<RT> S);

    template <class T> 
    void HermSV_Decompose(MatrixView<T> U, DiagMatrixView<RT> S);
    template <class T> 
    void SymSV_Decompose(
        MatrixView<T> U,
        DiagMatrixView<RT> S, MVP<T> Vt, RT& logdet, T& signdet);
    template <class T> 
    void SV_Decompose(
        SymMatrixView<T> A, DiagMatrixView<RT> S, MVP<T> Vt);

    template <class T, class T1> 
    void HermSV_Inverse(
        const GenMatrix<T1>& U, const GenDiagMatrix<RT1>& S, 
        int kmax, SymMatrixView<T> sinv);
    template <class T, class T1> 
    void SymSV_Inverse(
        const GenMatrix<T1>& U, const GenDiagMatrix<RT1>& S,
        const GenMatrix<T1>& Vt, int kmax, SymMatrixView<T> sinv);

    template <class T> 
    void HermTridiagonalChopSmallElements(
        VectorView<T> D, VectorView<T> E);

    template <class T> 
    void EigenFromTridiagonal_QR(
        MVP<T> U, VectorView<RT> D, VectorView<RT> E);

    template <class T> 
    void EigenFromTridiagonal_DC(
        MVP<T> U, VectorView<RT> D, VectorView<RT> E, bool UisI);

    template <class T> 
    void FindDCEigenValues(
        Vector<T>& S, const T rho, const GenVector<T>& D, const GenVector<T>& z,
        Matrix<T,ColMajor>& diffmat);

    template <class T> 
    void FindDCEigenValues(
        Vector<T>& S, const T rho, const GenVector<T>& D, const GenVector<T>& z);

#undef RT

#define CT std::complex<T>

    // Again, Microsoft Visual C++ needs this extra parameter to 
    // get the resolution of the overload right.
    template <class T, class T1> 
    inline void CallHermSV_Inverse(
        T , const GenMatrix<T1>& U, const GenDiagMatrix<RT1>& S,
        int kmax, SymMatrixView<T> sinv)
    { HermSV_Inverse(U,S,kmax,sinv); }

    template <class T, class T1> 
    inline void CallHermSV_Inverse(
        T , const GenMatrix<T1>& U, const GenDiagMatrix<RT1>& S, 
        int kmax, SymMatrixView<CT> sinv)
    { HermSV_Inverse(U,S,kmax,sinv); }

    template <class T, class T1> 
    inline void CallSymSV_Inverse(
        T , const GenMatrix<T1>& U, const GenDiagMatrix<RT1>& S,
        const GenMatrix<T1>& Vt, int kmax, SymMatrixView<T> sinv)
    { SymSV_Inverse(U,S,Vt,kmax,sinv); }

    template <class T, class T1> 
    inline void CallSymSV_Inverse(
        T , const GenMatrix<T1>& U, const GenDiagMatrix<RT1>& S,
        const GenMatrix<T1>& Vt, int kmax, SymMatrixView<CT> sinv)
    { SymSV_Inverse(U,S,Vt,kmax,sinv); }

    // Specialize disallowed complex combinations:
    template <class T>
    inline void CallHermSV_Inverse(
        CT , const GenMatrix<CT>& , const GenDiagMatrix<T>& , 
        int , SymMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T>
    inline void CallSymSV_Inverse(
        CT , const GenMatrix<CT>& , const GenDiagMatrix<T>& ,
        const GenMatrix<CT>& , int , SymMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

#undef CT
 
}

#endif
