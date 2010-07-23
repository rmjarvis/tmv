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
#ifndef TMV_SVDiv_H
#define TMV_SVDiv_H

#include "tmv/TMV_BaseMatrix.h"
#include "tmv/TMV_BaseDiagMatrix.h"

namespace tmv {

#define RT TMV_RealType(T)

    template <class T> 
    struct MVP 
    {
        MVP(const MatrixView<T>& m) : mp(&m) {}
        MVP(const MatrixView<T>* m) : mp(m) {}
        MVP(int) : mp(0) {}
        operator const MatrixView<T>*() { return mp; }
        const MatrixView<T>* operator->() { return mp; }
        const MatrixView<T>& operator*() { return *mp; }

        const MatrixView<T>* mp;
    };

    template <class T> 
    void SV_Decompose(const MatrixView<T>& U, const DiagMatrixView<RT>& S, 
                      MVP<T> V, RT& logdet, T& signdet, bool StoreU=false);

    template <class T, class Tm, class Tx> 
    void SV_LDiv(const GenMatrix<T>& U, const GenDiagMatrix<RT>& S,
                 const GenMatrix<T>& V, int kmax,
                 const GenMatrix<Tm>& m, const MatrixView<Tx>& x);
    template <class T, class Tm, class Tx> 
    void SV_RDiv(const GenMatrix<T>& U, const GenDiagMatrix<RT>& S,
                 const GenMatrix<T>& V, int kmax,
                 const GenMatrix<Tm>& m, const MatrixView<Tx>& x);


    template <class T> 
    void Bidiagonalize(const MatrixView<T>& A, const VectorView<T>& Ubeta,
                       const VectorView<T>& Vbeta, const VectorView<RT>& D,
                       const VectorView<RT>& E, T& signdet);

    template <class T> 
    void BidiagonalChopSmallElements(
        const VectorView<T>& D, const VectorView<T>& E, bool* zd=0);

    template <class T> 
    void BidiagonalZeroFirstRow(
        MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E);

    template <class T> 
    void BidiagonalZeroLastCol(
        const VectorView<RT>& D, const VectorView<RT>& E, MVP<T> V);

    template <class T> 
    void SV_DecomposeFromBidiagonal(
        MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E, MVP<T> V,
        bool SetUV=false);

    template <class T> 
    void DoSVDecomposeFromBidiagonal(
        MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E, MVP<T> V,
        bool UisI, bool VisI);

    template <class T> 
    void SV_DecomposeFromBidiagonal_QR(
        MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E, MVP<T> V);

    template <class T> 
    void SV_DecomposeFromBidiagonal_DC(
        MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E, MVP<T> V,
        bool UisI, bool VisI);

    template <class T> 
    void FindDCSingularValues(
        Vector<T>& S, const T rho, const GenVector<T>& D, const GenVector<T>& z,
        Matrix<T,ColMajor>& diff);

    template <class T> 
    void FindDCSingularValues(
        Vector<T>& S, const T rho, const GenVector<T>& D, const GenVector<T>& z);

#undef RT

} // namespace mv

#endif
