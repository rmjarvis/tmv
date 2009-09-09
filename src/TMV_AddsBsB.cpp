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



#include "TMV_Blas.h"
#include "TMV_SymBandMatrixArithFunc.h"
#include "TMV_SymBandMatrix.h"
#include "TMV_BandMatrixArith.h"
#include "TMV_TriMatrix.h"
#include "TMV_TriMatrixArith.h"
#include "TMV_SymBandMatrixArith.h"

//#define XDEBUG

#ifdef XDEBUG
#include "TMV_MatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {


  template <class T, class Ta> void AddMM(
      const T alpha, const GenSymBandMatrix<Ta>& A, const BandMatrixView<T>& B)
  {
#ifdef XTEST
    TMVAssert(A.HermOK());
#endif

#ifdef XDEBUG
    Matrix<Ta> A0 = A;
    Matrix<T> B0 = B;
    Matrix<T> B2 = alpha*A0 + B0;
    //cout<<"Start SymAddMM: alpha = "<<alpha<<endl;
    //cout<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cout<<"B = "<<Type(B)<<"  "<<B0<<endl;
#endif

    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == B.rowsize());
    TMVAssert(B.nlo() >= A.nlo());
    TMVAssert(B.nhi() >= A.nlo());
    if (A.size() > 0) {
      if (B.isconj()) AddMM(CONJ(alpha),A.Conjugate(),B.Conjugate());
      else {
	if (SameStorage(A,B)) {
	  if (B.isrm()) {
	    BandMatrix<Ta,RowMajor> tempA(A);
	    AddMM(alpha,tempA,B);
	  } else {
	    BandMatrix<Ta,ColMajor> tempA(A);
	    AddMM(alpha,tempA,B);
	  }
	}
	else {
	  AddMM(alpha,A.UpperBand(),BandMatrixView<T>(B,0,A.nlo()));
	  if (A.nlo()>0)
	    AddMM(alpha,A.LowerBandOff(),
		BandMatrixView<T>(B,A.nlo(),0).Diags(-A.nlo(),0));
	}
      }
    }
#ifdef XDEBUG
    //cout<<"Done\n";
    if (Norm(B-B2) > 0.001*(Norm(B0)+ABS(alpha)*Norm(A))) {
      cerr<<"SymAddMM: alpha = "<<alpha<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"->B = "<<B<<endl;
      cerr<<"B2 = "<<B2<<endl;
      abort();
    }
#endif
  }


  template <class T, class Ta, class Tb> void AddMM(
      const T alpha, const GenSymBandMatrix<Ta>& A,
      const T beta, const GenSymBandMatrix<Tb>& B, const BandMatrixView<T>& C)
  { 
#ifdef XTEST
    TMVAssert(A.HermOK());
    TMVAssert(B.HermOK());
#endif
#ifdef XDEBUG
    Matrix<Ta> A0 = A;
    Matrix<Tb> B0 = B;
    Matrix<T> C2 = alpha*A0 + beta*B0;
    //cout<<"Start SymAddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
    //cout<<"A = "<<Type(A)<<"  "<<A0<<endl;
    //cout<<"B = "<<Type(B)<<"  "<<B0<<endl;
#endif

    TMVAssert(A.size() == B.size());
    TMVAssert(C.rowsize() == A.size());
    TMVAssert(C.colsize() == A.size());

    if (A.size() > 0) {
      if (SameStorage(A,C)) {
	if (SameStorage(B,C)) {
	  if (C.isrm()) {
	    BandMatrix<T,RowMajor> tempC = beta*B;
	    tempC += alpha*A;
	    C = tempC;
	  } else {
	    BandMatrix<T,ColMajor> tempC = beta*B;
	    tempC += alpha*A;
	    C = tempC;
	  }
	} else {
	  C = alpha*A;
	  AddMM(beta,B,C);
	}
      } else {
	C = beta*B;
	AddMM(alpha,A,C);
      }
    }

#ifdef XDEBUG
    //cout<<"Done\n";
    if (Norm(C-C2) > 0.001*(ABS(alpha)*Norm(A0)+ABS(beta)*Norm(B0))) {
      cerr<<"Start SymAddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"->C = "<<Type(C)<<"  "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  } 

  template <class T, class Ta, class Tb> void AddMM(
      const T alpha, const GenSymBandMatrix<Ta>& A,
      const T beta, const GenSymBandMatrix<Tb>& B, const MatrixView<T>& C)
  { 
#ifdef XTEST
    TMVAssert(A.HermOK());
    TMVAssert(B.HermOK());
#endif
#ifdef XDEBUG
    Matrix<Ta> A0 = A;
    Matrix<Tb> B0 = B;
    Matrix<T> C2 = alpha*A0 + beta*B0;
    //cout<<"Start SymAddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
    //cout<<"A = "<<Type(A)<<"  "<<A0<<endl;
    //cout<<"B = "<<Type(B)<<"  "<<B0<<endl;
#endif

    TMVAssert(A.size() == B.size());
    TMVAssert(C.rowsize() == A.size());
    TMVAssert(C.colsize() == A.size());

    const int N = A.size();
    const int k = MAX(A.nlo(),B.nlo());

    if (N > 0) {
      if (SameStorage(A,C) || SameStorage(B,C)) {
	AddMM(alpha,A,beta,B,BandMatrixView<T>(C,k,k));
	UpperTriMatrixViewOf(C.SubMatrix(0,N-k-1,k+1,N)).Zero();
	LowerTriMatrixViewOf(C.SubMatrix(k+1,N,0,N-k-1)).Zero();
      } else {
	C.Zero();
	AddMM(alpha,A,beta,B,BandMatrixView<T>(C,k,k));
      }
    }

#ifdef XDEBUG
    //cout<<"Done\n";
    if (Norm(C-C2) > 0.001*(ABS(alpha)*Norm(A0)+ABS(beta)*Norm(B0))) {
      cerr<<"Start SymAddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"->C = "<<Type(C)<<"  "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta, class Tb> void AddMM(
      const T alpha, const GenSymBandMatrix<Ta>& A,
      const T beta, const GenBandMatrix<Tb>& B, const BandMatrixView<T>& C)
  { 
#ifdef XTEST
    TMVAssert(A.HermOK());
#endif
#ifdef XDEBUG
    Matrix<Ta> A0 = A;
    Matrix<Tb> B0 = B;
    Matrix<T> C2 = alpha*A0 + beta*B0;
    //cout<<"Start SymAddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
    //cout<<"A = "<<Type(A)<<"  "<<A0<<endl;
    //cout<<"B = "<<Type(B)<<"  "<<B0<<endl;
#endif

    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == B.rowsize());
    TMVAssert(C.rowsize() == A.size());
    TMVAssert(C.colsize() == A.size());
    TMVAssert(C.nlo() >= A.nlo());
    TMVAssert(C.nlo() >= B.nlo());
    TMVAssert(C.nhi() >= A.nhi());
    TMVAssert(C.nhi() >= B.nhi());

    if (A.size() > 0) {
      if (SameStorage(A,C)) {
	if (SameStorage(B,C)) {
	  if (C.isrm()) {
	    BandMatrix<T,RowMajor> tempC(C.colsize(),C.rowsize(),
		C.nlo(),C.nhi());
	    tempC = beta*B;
	    tempC += alpha*A;
	    C = tempC;
	  } else {
	    BandMatrix<T,ColMajor> tempC(C.colsize(),C.rowsize(),
		C.nlo(),C.nhi());
	    tempC = beta*B;
	    tempC += alpha*A;
	    C = tempC;
	  }
	} else {
	  C = alpha*A;
	  C += beta*B;
	}
      } else {
	C = beta*B;
	AddMM(alpha,A,C);
      }
    }

#ifdef XDEBUG
    //cout<<"Done\n";
    if (Norm(C-C2) > 0.001*(ABS(alpha)*Norm(A0)+ABS(beta)*Norm(B0))) {
      cerr<<"SymAddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"->C = "<<Type(C)<<"  "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  } 

  template <class T, class Ta, class Tb> void AddMM(
      const T alpha, const GenSymBandMatrix<Ta>& A,
      const T beta, const GenMatrix<Tb>& B, const MatrixView<T>& C)
  { 
#ifdef XTEST
    TMVAssert(A.HermOK());
#endif
#ifdef XDEBUG
    Matrix<Ta> A0 = A;
    Matrix<Tb> B0 = B;
    Matrix<T> C2 = alpha*A0 + beta*B0;
    //cout<<"Start SymAddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
    //cout<<"A = "<<Type(A)<<"  "<<A0<<endl;
    //cout<<"B = "<<Type(B)<<"  "<<B0<<endl;
#endif

    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == B.rowsize());
    TMVAssert(C.rowsize() == A.size());
    TMVAssert(C.colsize() == A.size());

    const int N = A.size();
    const int k = A.nlo();
    if (N > 0) {
      AddMM(alpha,A,beta,ConstBandMatrixView<Tb>(B,k,k),
	 BandMatrixView<T>(C,k,k));
      UpperTriMatrixViewOf(C.SubMatrix(0,N-k-1,k+1,N)) =
	beta * UpperTriMatrixViewOf(B.SubMatrix(0,N-k-1,k+1,N)); 
      LowerTriMatrixViewOf(C.SubMatrix(k+1,N,0,N-k-1)) =
	beta * LowerTriMatrixViewOf(B.SubMatrix(k+1,N,0,N-k-1));
    }

#ifdef XDEBUG
    //cout<<"Done\n";
    if (Norm(C-C2) > 0.001*(ABS(alpha)*Norm(A0)+ABS(beta)*Norm(B0))) {
      cerr<<"SymAddMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"->C = "<<Type(C)<<"  "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_AddsBsB.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


