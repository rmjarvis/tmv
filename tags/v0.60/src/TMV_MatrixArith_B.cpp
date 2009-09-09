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
#include "TMV_VectorArith.h"
#include "TMV_MatrixArith.h"

//#define XDEBUG

#ifdef XDEBUG
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

  //
  // MultXM
  //

  template <class T, class Ta> inline void RowMajorMultXM(
      const Ta alpha, const MatrixView<T>& A)
  {
    TMVAssert(A.isrm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.colsize() > 0 && A.rowsize() > 0);
    TMVAssert(alpha != Ta(1));
    TMVAssert(alpha != Ta(0));

    const size_t M = A.colsize();
    const size_t N = A.rowsize();

    T* Ai0 = A.ptr();

    for(size_t i=M;i>0;--i,Ai0+=A.stepi()) {
      // A.row(i) *= alpha;
      T* Aij = Ai0;
      for(size_t j=N;j>0;--j,++Aij) *Aij *= alpha;
    }
  }

  template <class T> void MultXM(const T alpha, const MatrixView<T>& A)
    // A = alpha * A
  {
#ifdef XDEBUG
    Matrix<T> A0 = A;
    Matrix<T> A2 = A;
    for(size_t i=0;i<A.colsize();i++)
      for(size_t j=0;j<A.rowsize();j++)
	A2(i,j) *= alpha;
    //cerr<<"MultXM: alpha = "<<alpha<<", A = "<<Type(A)<<"  "<<A<<endl;
#endif

    if (A.colsize() > 0 && A.rowsize() > 0 && alpha != T(1)) {
      if (A.isconj()) MultXM(CONJ(alpha),A.Conjugate());
      else if (alpha == T(0)) A.Zero();
      else if (A.CanLinearize()) A.LinearView() *= alpha;
      else if (A.isrm())
	if (IsComplex(T()) && IMAG(alpha) == RealType(T)(0))
	  RowMajorMultXM(REAL(alpha),A);
	else
	  RowMajorMultXM(alpha,A);
      else if (A.iscm())
	if (IsComplex(T()) && IMAG(alpha) == RealType(T)(0))
	  RowMajorMultXM(REAL(alpha),A.Transpose());
	else 
	  RowMajorMultXM(alpha,A.Transpose());
      else 
	if (A.colsize() < A.rowsize())
	  for(size_t i=0;i<A.colsize();++i) A.row(i) *= alpha;
	else 
	  for(size_t j=0;j<A.colsize();++j) A.col(j) *= alpha;
    }
#ifdef XDEBUG
    //cerr<<"Done: A = "<<A<<endl;
    // Doing A-A2 becomes recursive call to MultXM
    Matrix<T> diff(A.colsize(),A.rowsize());
    for(size_t i=0;i<A.colsize();i++)
      for(size_t j=0;j<A.rowsize();j++)
	diff(i,j) = A(i,j) - A2(i,j);
    if (Norm(diff) > 0.001*ABS(alpha)*Norm(A)) {
      cerr<<"TriMultXM: alpha = "<<alpha<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> A = "<<A<<endl;
      cerr<<"A2 = "<<A2<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta> void ElementProd(
      const T alpha, const GenMatrix<Ta>& A, const MatrixView<T>& B)
  {
    TMVAssert(A.colsize() == B.colsize());
    TMVAssert(A.rowsize() == B.rowsize());
    if (A.stor() == B.stor() && A.CanLinearize() && B.CanLinearize()) {
      TMVAssert(A.stepi() == B.stepi() && A.stepj() == B.stepj());
      ElementProd(alpha,A.ConstLinearView(),B.LinearView());
    } else if (B.isrm()) {
      for(size_t i=0;i<B.colsize();i++)
	ElementProd(alpha,A.row(i),B.row(i));
    } else {
      for(size_t j=0;j<B.rowsize();j++)
	ElementProd(alpha,A.col(j),B.col(j));
    }
  }

  template <class T, class Ta, class Tb> void AddElementProd(
      const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == C.rowsize());
    TMVAssert(B.colsize() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    if (A.stor() == C.stor() && B.stor() == C.stor() &&
	A.CanLinearize() && B.CanLinearize() && C.CanLinearize()) {
      TMVAssert(A.stepi() == C.stepi() && A.stepj() == C.stepj());
      TMVAssert(B.stepi() == C.stepi() && B.stepj() == C.stepj());
      AddElementProd(alpha,A.ConstLinearView(),B.ConstLinearView(),
	  C.LinearView());
    } else if (C.isrm()) {
      for(size_t i=0;i<C.colsize();i++)
	AddElementProd(alpha,A.row(i),B.row(i),C.row(i));
    } else {
      for(size_t j=0;j<C.rowsize();j++)
	AddElementProd(alpha,A.col(j),B.col(j),C.col(j));
    }
  }

#define InstFile "TMV_MatrixArith_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


