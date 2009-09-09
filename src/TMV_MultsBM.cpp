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



#include "TMV_SymBandMatrixArithFunc.h"
#include "TMV_SymBandMatrix.h"
#include "TMV_Matrix.h"
#include "TMV_BandMatrixArith.h"

//#define XDEBUG

#ifdef XDEBUG
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define SYM_MM_BLOCKSIZE TMV_BLOCKSIZE
#define SYM_MM_BLOCKSIZE2 (TMV_BLOCKSIZE/2)
#else
#define SYM_MM_BLOCKSIZE 64
#define SYM_MM_BLOCKSIZE2 32
#endif

  //
  // MultMM
  //

  template <bool add, class T, class Ta, class Tb> static void DoMultMM(
      const T alpha, const GenSymBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));

    if (add) C += alpha * A.LowerBand() * B;
    else C = alpha * A.LowerBand() * B;

    const int N = A.size();
    if (N > 1 && A.nlo() > 0)
      C.Rows(0,N-1) += A.UpperBandOff() * B.Rows(1,N);
  }

  template <bool add, class T, class Ta, class Tb> static void FullTempMultMM(
      const T alpha, const GenSymBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
    if (C.isrm()) {
      Matrix<T,RowMajor> C2(C.colsize(),C.rowsize());
      DoMultMM<false>(T(1),A,B,C2.View());
      if (add) C += alpha*C2;
      else C = alpha*C2;
    } else {
      Matrix<T,ColMajor> C2(C.colsize(),C.rowsize());
      DoMultMM<false>(T(1),A,B,C2.View());
      if (add) C += alpha*C2;
      else C = alpha*C2;
    }
  }

  template <bool add, class T, class Ta, class Tb> static void BlockTempMultMM(
      const T alpha, const GenSymBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
    const int N = C.rowsize();
    for(int j=0;j<N;) {
      int j2 = MIN(N,j+SYM_MM_BLOCKSIZE);
      if (IMAG(alpha) == RealType(T)(0)) {
	if (C.isrm()) {
	  Matrix<Tb,RowMajor> B2 = REAL(alpha) * B.Cols(j,j2);
	  DoMultMM<add>(T(1),A,B2,C.Cols(j,j2));
	} else {
	  Matrix<Tb,ColMajor> B2 = REAL(alpha) * B.Cols(j,j2);
	  DoMultMM<add>(T(1),A,B2,C.Cols(j,j2));
	}
      } else {
	if (C.isrm()) {
	  Matrix<T,RowMajor> B2 = alpha * B.Cols(j,j2);
	  DoMultMM<add>(T(1),A,B2,C.Cols(j,j2));
	} else {
	  Matrix<T,ColMajor> B2 = alpha * B.Cols(j,j2);
	  DoMultMM<add>(T(1),A,B2,C.Cols(j,j2));
	}
      }
      j = j2;
    }
  }

  template <bool add, class T, class Ta, class Tb> void MultMM(
      const T alpha, const GenSymBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
    // C (+)= alpha * A * B
  {
#ifdef XTEST
    TMVAssert(A.HermOK());
#endif
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
#ifdef XDEBUG
    //cout<<"Start MultMM: alpha = "<<alpha<<endl;
    //cout<<"A = "<<A.cptr()<<"  "<<Type(A)<<"  "<<A<<endl;
    //cout<<"B = "<<B.cptr()<<"  "<<Type(B)<<"  "<<B<<endl;
    //cout<<"C = "<<C.cptr()<<"  "<<Type(C)<<"  "<<C<<endl;
    Matrix<Ta> A0 = A;
    Matrix<Tb> B0 = B;
    Matrix<T> C0 = C;
    Matrix<T> C2 = alpha*A0*B0;
    if (add) C2 += C0;
#endif

    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (alpha == T(0)) {
	if (!add) C.Zero();
      }
      else if (C.isconj()) {
	MultMM<add>(CONJ(alpha),A.Conjugate(),B.Conjugate(),C.Conjugate());
      }
      else if (SameStorage(A,C)) 
	FullTempMultMM<add>(alpha,A,B,C);
      else if (SameStorage(B,C)) 
	if (C.stepi() == B.stepi() && C.stepj() == B.stepj())
	  BlockTempMultMM<add>(alpha,A,B,C);
	else
	  FullTempMultMM<add>(alpha,A,B,C);
      else BlockTempMultMM<add>(alpha, A, B, C);
    }

#ifdef XDEBUG
    //cout<<"Done: C = "<<C<<endl;
    if (Norm(C-C2) > 0.001*(ABS(alpha)*Norm(A0)*Norm(B0)+
	  (add?Norm(C0):RealType(T)(0)))) {
      cerr<<"MultMM: alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"--> C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_MultsBM.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


