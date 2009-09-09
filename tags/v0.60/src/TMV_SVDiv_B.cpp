///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2007                                                        //
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
// along with this program in the file gpl.txt.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



#include "TMV_Matrix.h"
#include "TMV_DiagMatrix.h"
#include "TMV_SVDiv.h"
#include "TMV_MatrixArith.h"
#include "TMV_DiagMatrixArith.h"

//#define XDEBUG

#ifdef XDEBUG
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

  //
  // LDiv
  //

  template <class T, class T1, class T2> void SV_LDiv(
      const GenMatrix<T1>& U, const GenVector<RealType(T1)>& S, 
      const GenMatrix<T1>& V, size_t kmax,
      const GenMatrix<T2>& m, const MatrixView<T>& x)
  {
    // A x = m
    // U S V x = m
    // x = Vt S^-1 Ut m
    TMVAssert(m.colsize() == U.colsize()); // = M
    TMVAssert(x.colsize() == V.rowsize()); // = N
    TMVAssert(x.rowsize() == m.rowsize()); // = R
    TMVAssert(kmax <= V.rowsize()); // = K
    TMVAssert(kmax <= U.colsize());
#ifdef XDEBUG
    Matrix<T1> A = U * DiagMatrixViewOf(S) * V;
    Matrix<T2> m0(m);
#endif

    Matrix<T> m2 = U.Adjoint().Rows(0,kmax) * m; // KxR
    m2 /= DiagMatrixViewOf(S.SubVector(0,kmax));
    x = V.Adjoint().Cols(0,kmax) * m2; // NxR

#ifdef XDEBUG
    // Note: this test only works for square matrices
    Matrix<T> mm = A*x;
    if (U.IsSquare() && Norm(m0-mm) > 0.001 * Norm(A)*Norm(m0)) {
      cerr<<"SV_LDiv\n";
      cerr<<"U = "<<U<<endl;
      cerr<<"S = "<<S<<endl;
      cerr<<"V = "<<V<<endl;
      cerr<<"A = USV = "<<A<<endl;
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

  template <class T, class T1, class T2> void SV_RDiv(
      const GenMatrix<T1>& U, const GenVector<RealType(T1)>& S, 
      const GenMatrix<T1>& V, size_t kmax,
      const GenMatrix<T2>& m, const MatrixView<T>& x) 
  {
    // x A = m
    // x U S V = m
    // x = m Vt S^-1 Ut
    TMVAssert(m.rowsize() == V.rowsize()); // = N
    TMVAssert(x.rowsize() == U.colsize()); // = M
    TMVAssert(x.colsize() == m.colsize()); // = R
    TMVAssert(kmax <= U.colsize()); // = K
    TMVAssert(kmax <= V.rowsize());
#ifdef XDEBUG
    Matrix<T1> A = U * DiagMatrixViewOf(S) * V;
    Matrix<T2> m0(m);
#endif

    Matrix<T,ColMajor> m2 = m * V.Adjoint().Cols(0,kmax); // = RxK
    m2 %= DiagMatrixViewOf(S.SubVector(0,kmax));
    x = m2 * U.Adjoint().Rows(0,kmax); // = RxM

#ifdef XDEBUG
    // Note: this test only works for square matrices
    Matrix<T> mm = x*A;
    if (U.IsSquare() && Norm(m0-mm) > 0.001 * Norm(A)*Norm(m0)) {
      cerr<<"SV_RDiv\n";
      cerr<<"U = "<<U<<endl;
      cerr<<"S = "<<S<<endl;
      cerr<<"V = "<<V<<endl;
      cerr<<"A = USV = "<<A<<endl;
      cerr<<"m0 = "<<m0<<endl;
      cerr<<"x = "<<x<<endl;
      cerr<<"xA = "<<mm<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_SVDiv_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


