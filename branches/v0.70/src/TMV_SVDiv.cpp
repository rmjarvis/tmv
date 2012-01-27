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


//#define XDEBUG


#include "TMV_SVDiv.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_DiagMatrixArith.h"
#include "tmv/TMV_MatrixArith.h"

#ifdef XDEBUG
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

#define RT TMV_RealType(T)

    //
    // LDiv
    //

    template <class T, class Tm, class Tx> 
    void SV_LDiv(
        const GenMatrix<T>& U, const GenDiagMatrix<RT>& SS, 
        const GenMatrix<T>& Vt, int kmax,
        const GenMatrix<Tm>& m, MatrixView<Tx> x)
    {
        // A x = m
        // U S Vt x = m
        // x = V S^-1 Ut m
        TMVAssert(m.colsize() == U.colsize()); // = M
        TMVAssert(x.colsize() == Vt.rowsize()); // = N
        TMVAssert(x.rowsize() == m.rowsize()); // = R
        TMVAssert(kmax <= Vt.rowsize()); // = K
        TMVAssert(kmax <= U.colsize());
#ifdef XDEBUG
        Matrix<T> A = U * SS * Vt;
        Matrix<Tm> m0(m);
#endif

        Matrix<Tx,RowMajor> m2 = U.adjoint().rowRange(0,kmax) * m; // KxR
        m2 /= SS.subDiagMatrix(0,kmax);
        x = Vt.adjoint().colRange(0,kmax) * m2; // NxR

#ifdef XDEBUG
        // Note: this test only works for square matrices
        Matrix<Tx> mm = A*x;
        if (U.isSquare() && kmax==SS.size() && 
            !(Norm(m0-mm) <= 0.001 * Norm(A)*Norm(m0))) {
            cerr<<"SV_LDiv\n";
            cerr<<"U = "<<U<<endl;
            cerr<<"S = "<<SS.diag()<<endl;
            cerr<<"Vt = "<<Vt<<endl;
            cerr<<"A = USVt = "<<A<<endl;
            cerr<<"m0 = "<<m0<<endl;
            cerr<<"x = "<<x<<endl;
            cerr<<"Ax = "<<mm<<endl;
            abort();
        }
#endif
    }

    //
    // RDiv
    //

    template <class T, class Tm, class Tx> 
    void SV_RDiv(
        const GenMatrix<T>& U, const GenDiagMatrix<RT>& SS, 
        const GenMatrix<T>& Vt, int kmax,
        const GenMatrix<Tm>& m, MatrixView<Tx> x) 
    {
        // x A = m
        // x U S Vt = m
        // x = m V S^-1 Ut
        TMVAssert(m.rowsize() == Vt.rowsize()); // = N
        TMVAssert(x.rowsize() == U.colsize()); // = M
        TMVAssert(x.colsize() == m.colsize()); // = R
        TMVAssert(kmax <= U.colsize()); // = K
        TMVAssert(kmax <= Vt.rowsize());
#ifdef XDEBUG
        Matrix<T> A = U * SS * Vt;
        Matrix<Tm> m0(m);
#endif

        Matrix<Tx,RowMajor> m2 = m * Vt.adjoint().colRange(0,kmax); // = RxK
        m2 %= SS.subDiagMatrix(0,kmax);
        x = m2 * U.adjoint().rowRange(0,kmax); // = RxM

#ifdef XDEBUG
        // Note: this test only works for square matrices
        Matrix<Tx> mm = x*A;
        if (U.isSquare() && kmax==SS.size() && 
            !(Norm(m0-mm) > 0.001 * Norm(A)*Norm(m0))) {
            cerr<<"SV_RDiv\n";
            cerr<<"U = "<<U<<endl;
            cerr<<"S = "<<SS.diag()<<endl;
            cerr<<"Vt = "<<Vt<<endl;
            cerr<<"A = USVt = "<<A<<endl;
            cerr<<"m0 = "<<m0<<endl;
            cerr<<"x = "<<x<<endl;
            cerr<<"xA = "<<mm<<endl;
            abort();
        }
#endif
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_SVDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


