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
        MVP(MatrixView<T> m) : mp(&m) {}
        MVP(MatrixView<T>* m) : mp(m) {}
        MVP(int) : mp(0) {}
        operator MatrixView<T>*() { return mp; }
        MatrixView<T>* operator->() { return mp; }
        MatrixView<T> operator*() { return *mp; }

        MatrixView<T>* mp;
    };

    template <class T> 
    void SV_Decompose(
        MatrixView<T> U, DiagMatrixView<RT> S, 
        MVP<T> Vt, RT& logdet, T& signdet, bool StoreU=false);

    template <class T, class Tm, class Tx> 
    void SV_LDiv(
        const GenMatrix<T>& U, const GenDiagMatrix<RT>& S,
        const GenMatrix<T>& Vt, int kmax,
        const GenMatrix<Tm>& m, MatrixView<Tx> x);
    template <class T, class Tm, class Tx> 
    void SV_RDiv(
        const GenMatrix<T>& U, const GenDiagMatrix<RT>& S,
        const GenMatrix<T>& Vt, int kmax,
        const GenMatrix<Tm>& m, MatrixView<Tx> x);


    template <class T> 
    void Bidiagonalize(
        MatrixView<T> A, VectorView<T> Ubeta,
        VectorView<T> Vtbeta, VectorView<RT> D, VectorView<RT> E, T& signdet);

    template <class T> 
    void BidiagonalChopSmallElements(
        VectorView<T> D, VectorView<T> E, bool* zd=0);

    template <class T> 
    void BidiagonalZeroFirstRow(
        MVP<T> U, VectorView<RT> D, VectorView<RT> E);

    template <class T> 
    void BidiagonalZeroLastCol(
        VectorView<RT> D, VectorView<RT> E, MVP<T> Vt);

    template <class T> 
    void SV_DecomposeFromBidiagonal(
        MVP<T> U, VectorView<RT> D, VectorView<RT> E, MVP<T> Vt,
        bool SetUV=false);

    template <class T> 
    void DoSVDecomposeFromBidiagonal(
        MVP<T> U, VectorView<RT> D, VectorView<RT> E, MVP<T> Vt,
        bool UisI, bool VisI);

    template <class T> 
    void SV_DecomposeFromBidiagonal_QR(
        MVP<T> U, VectorView<RT> D, VectorView<RT> E, MVP<T> Vt);

    template <class T> 
    void SV_DecomposeFromBidiagonal_DC(
        MVP<T> U, VectorView<RT> D, VectorView<RT> E, MVP<T> Vt,
        bool UisI, bool VisI);

    template <class T> 
    void FindDCSingularValues(
        Vector<T>& S, const T rho, const GenVector<T>& D, const GenVector<T>& z,
        Matrix<T,ColMajor>& diff);

    template <class T> 
    void FindDCSingularValues(
        Vector<T>& S, const T rho, const GenVector<T>& D, const GenVector<T>& z);


#define CTx std::complex<Tx>

    // Some compilers (specifically MS Visual C++ at least) need this
    // dummy variable in the front to be able to resolve the overloads
    // between the allowed versions and the below disallowed versions.
    template <class T, class Tm, class Tx>
    inline void CallSV_LDiv(
        Tx , const GenMatrix<T>& U, const GenDiagMatrix<RT>& S,
        const GenMatrix<T>& Vt, int kmax,
        const GenMatrix<Tm>& m, MatrixView<Tx> x)
    { SV_LDiv(U,S,Vt,kmax,m,x); }

    template <class T, class Tm, class Tx>
    inline void CallSV_LDiv(
        Tx , const GenMatrix<T>& U, const GenDiagMatrix<RT>& S,
        const GenMatrix<T>& Vt, int kmax,
        const GenMatrix<Tm>& m, MatrixView<CTx> x)
    { SV_LDiv(U,S,Vt,kmax,m,x); }
   
    template <class T, class Tm, class Tx>
    inline void CallSV_RDiv(
        Tx , const GenMatrix<T>& U, const GenDiagMatrix<RT>& S,
        const GenMatrix<T>& Vt, int kmax,
        const GenMatrix<Tm>& m, MatrixView<Tx> x)
    { SV_RDiv(U,S,Vt,kmax,m,x); }
   
    template <class T, class Tm, class Tx>
    inline void CallSV_RDiv(
        Tx , const GenMatrix<T>& U, const GenDiagMatrix<RT>& S,
        const GenMatrix<T>& Vt, int kmax,
        const GenMatrix<Tm>& m, MatrixView<CTx> x)
    { SV_RDiv(U,S,Vt,kmax,m,x); }
   
#undef CTx
#undef RT

    // Specialize disallowed complex combinations:
#define CT std::complex<T>

    template <class T>
    inline void CallSV_LDiv(
        CT , const GenMatrix<CT>& , const GenDiagMatrix<T>& ,
        const GenMatrix<CT>& , int ,
        const GenMatrix<CT>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T>
    inline void CallSV_LDiv(
        CT , const GenMatrix<CT>& , const GenDiagMatrix<T>& ,
        const GenMatrix<CT>& , int ,
        const GenMatrix<T>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T>
    inline void CallSV_LDiv(
        CT , const GenMatrix<T>& , const GenDiagMatrix<T>& ,
        const GenMatrix<T>& , int ,
        const GenMatrix<CT>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T>
    inline void CallSV_RDiv(
        CT , const GenMatrix<CT>& , const GenDiagMatrix<T>& ,
        const GenMatrix<CT>& , int ,
        const GenMatrix<CT>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T>
    inline void CallSV_RDiv(
        CT , const GenMatrix<CT>& , const GenDiagMatrix<T>& ,
        const GenMatrix<CT>& , int ,
        const GenMatrix<T>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T>
    inline void CallSV_RDiv(
        CT , const GenMatrix<T>& , const GenDiagMatrix<T>& ,
        const GenMatrix<T>& , int ,
        const GenMatrix<CT>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

#undef CT

} // namespace mv

#endif
