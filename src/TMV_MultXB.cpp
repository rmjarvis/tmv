///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2008                                                        //
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



#include "TMV_BandMatrixArithFunc.h"
#include "TMV_BandMatrix.h"
#include "TMV_VectorArith.h"

//#define XDEBUG

#ifdef XDEBUG
#include "TMV_MatrixArith.h"
#include "TMV_BandMatrixArith.h"
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

  //
  // MultXM
  //

  template <class T> void MultXM(const T alpha, const BandMatrixView<T>& A)
    // A = alpha * A
  {
#ifdef XDEBUG
    Matrix<T> A0 = A;
    Matrix<T> A2 = alpha*A0;
#endif
    if (A.rowsize() > 0 && A.colsize() > 0 && alpha != T(1)) {
      if (A.isconj()) MultXM(CONJ(alpha),A.Conjugate());
      else if (alpha == T(0)) A.Zero();
      else if (A.CanLinearize()) A.LinearView() *= alpha;
      else 
	for(int i=-A.nlo();i<=A.nhi();++i) A.diag(i) *= alpha;
    }
#ifdef XDEBUG
    if (Norm(A2-A) > 0.001*ABS(alpha)*Norm(A0)) {
      cerr<<"MultXM: alpha = "<<alpha;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"A2 = "<<A2<<endl;
      cerr<<"Norm(diff) = "<<Norm(A-A2)<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta> void ElementProd(
      const T alpha, const GenBandMatrix<Ta>& A, const BandMatrixView<T>& B)
  {
    TMVAssert(A.colsize() == B.colsize());
    TMVAssert(A.rowsize() == B.rowsize());
    TMVAssert(A.nlo() == B.nlo());
    TMVAssert(A.nhi() == B.nhi());
    if (A.stor() == B.stor() && A.CanLinearize() && B.CanLinearize()) {
      TMVAssert(A.stepi() == B.stepi() && A.stepj() == B.stepj());
      ElementProd(alpha,A.ConstLinearView(),B.LinearView());
    } else {
      for(int i=-A.nlo();i<=A.nhi();++i) 
	ElementProd(alpha,A.diag(i),B.diag(i));
    }
  }

  template <class T, class Ta, class Tb> void AddElementProd(
      const T alpha, const GenBandMatrix<Ta>& A, const GenBandMatrix<Tb>& B,
      const BandMatrixView<T>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == C.rowsize());
    TMVAssert(B.colsize() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(A.nlo() == C.nlo());
    TMVAssert(A.nhi() == C.nhi());
    TMVAssert(B.nlo() == C.nlo());
    TMVAssert(B.nhi() == C.nhi());
    if (A.stor() == C.stor() && B.stor() == C.stor() &&
	A.CanLinearize() && B.CanLinearize() && C.CanLinearize()) {
      TMVAssert(A.stepi() == C.stepi() && A.stepj() == C.stepj());
      TMVAssert(B.stepi() == C.stepi() && B.stepj() == C.stepj());
      AddElementProd(alpha,A.ConstLinearView(),B.ConstLinearView(),
	  C.LinearView());
    } else {
      for(int i=-A.nlo();i<=A.nhi();++i) 
	AddElementProd(alpha,A.diag(i),B.diag(i),C.diag(i));
    }
  }

#define InstFile "TMV_MultXB.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


