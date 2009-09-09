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



#include "TMV_TriMatrix.h"
#include "TMV_VectorArith.h"
#include "TMV_TriMatrixArith.h"

//#define XDEBUG

#ifdef XDEBUG
#include "TMV_MatrixArith.h"
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

  //
  // MultXM
  //

  template <class T, class T1> inline void RowMajorMultXM(
      const T1 alpha, const UpperTriMatrixView<T>& A)
  {
    TMVAssert(A.isrm());
    TMVAssert(!A.isunit());
    TMVAssert(alpha != T1(1));
    TMVAssert(A.size() > 0);
    TMVAssert(A.ct() == NonConj);

    T* Aii = A.ptr();
    const int ds = A.stepi()+1;
    const size_t N = A.size();

    for(size_t len=N;len>0;--len,Aii+=ds) {
      // A.row(i,i,N) *= alpha;
      T* Aij = Aii;
      for(size_t j=len;j>0;--j,++Aij) *Aij *= alpha;
    }
  }

  template <class T, class T1> inline void ColMajorMultXM(
      const T1 alpha, const UpperTriMatrixView<T>& A)
  {
    TMVAssert(A.iscm());
    TMVAssert(!A.isunit());
    TMVAssert(alpha != T1(1));
    TMVAssert(A.size() > 0);
    TMVAssert(A.ct() == NonConj);

    T* A0j = A.ptr();
    const int Astepj = A.stepj();
    const size_t N = A.size();

    for(size_t j=N,len=1;j>0;--j,++len,A0j+=Astepj) {
      // A.col(j,0,j+1) *= alpha;
      T* Aij = A0j;
      for(size_t i=len;i>0;--i,++Aij) *Aij *= alpha;
    }
  }

  template <class T> void MultXM(const T alpha, const UpperTriMatrixView<T>& A)
    // A = alpha * A
  {
#ifdef XDEBUG
    Matrix<T> A0 = A;
    Matrix<T> A2 = alpha * A0;
    //cerr<<"MultXM: alpha = "<<alpha<<", A = "<<Type(A)<<" "<<A<<endl;
#endif

    if (A.size() > 0 && alpha != T(1)) {
      TMVAssert(!A.isunit());
      if (A.isconj()) MultXM(CONJ(alpha),A.Conjugate());
      else if (alpha == T(0)) A.Zero();
      else if (A.isrm())
	if (IMAG(alpha) == RealType(T)(0))
	  RowMajorMultXM(REAL(alpha),A);
	else
	  RowMajorMultXM(alpha,A);
      else if (A.iscm())
	if (IMAG(alpha) == RealType(T)(0))
	  ColMajorMultXM(REAL(alpha),A);
	else
	  ColMajorMultXM(alpha,A);
      else 
	for(size_t i=0;i<A.colsize();++i) 
	  A.row(i,i,A.rowsize()) *= alpha;
    }
#ifdef XDEBUG
    //cerr<<"Done MultXM: A = "<<A<<endl;
    if (Norm(Matrix<T>(A)-A2) > 0.001*ABS(alpha)*Norm(A)) {
      cerr<<"TriMultXM: alpha = "<<alpha<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> A = "<<A<<endl;
      cerr<<"A2 = "<<A2<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta> void ElementProd(
      const T alpha, const GenUpperTriMatrix<Ta>& A,
      const UpperTriMatrixView<T>& B)
  {
    TMVAssert(A.size() == B.size());
    const size_t N = B.size();
    if (B.isrm()) {
      for(size_t i=0;i<B.colsize();i++)
	ElementProd(alpha,A.row(i,i,N),B.row(i,i,N));
    } else {
      for(size_t j=0;j<B.rowsize();j++)
	ElementProd(alpha,A.col(j,0,j+1),B.col(j,0,j+1));
    }
  }

  template <class T, class Ta, class Tb> void AddElementProd(
      const T alpha, const GenUpperTriMatrix<Ta>& A,
      const GenUpperTriMatrix<Tb>& B, const UpperTriMatrixView<T>& C)
  {
    TMVAssert(A.size() == C.size());
    TMVAssert(B.size() == C.size());
    const size_t N = C.size();
    if (C.isrm()) {
      for(size_t i=0;i<C.colsize();i++)
	AddElementProd(alpha,A.row(i,i,N),B.row(i,i,N),C.row(i,i,N));
    } else {
      for(size_t j=0;j<C.rowsize();j++)
	AddElementProd(alpha,A.col(j,0,j+1),B.col(j,0,j+1),C.col(j,0,j+1));
    }
  }

#define InstFile "TMV_TriMatrixArith_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv

