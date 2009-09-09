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
#include "TMV_MatrixArith.h"
#include "TMV_TriMatrixArith.h"

//#define XDEBUG

#ifdef XDEBUG
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define TRI_MM_BLOCKSIZE TMV_BLOCKSIZE
#define TRI_MM_BLOCKSIZE2 (TMV_BLOCKSIZE/2)
#else
#define TRI_MM_BLOCKSIZE 64
#define TRI_MM_BLOCKSIZE2 32
#endif

  //
  // MultMM: M = U * L
  //

  template <bool add, class T, class Ta, class Tb> inline void ColMultMM(
      const T alpha, const GenUpperTriMatrix<Ta>& A,
      const GenLowerTriMatrix<Tb>& B, const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"MultMM: "<<alpha<<"  "<<Type(A)<<"  "<<A<<"  "<<Type(B)<<"  "<<B<<endl;
    //cerr<<Type(C)<<"  "<<C<<endl;
    Matrix<T> C2 = C;
    Matrix<T> C0 = C;
    Matrix<Ta> A0 = A;
    Matrix<Tb> B0 = B;
    if (add) C2 += alpha*A0*B0;
    else C2 = alpha*A0*B0;
    //cerr<<"ColMultMM UL\n";
    //cerr<<"A = "<<A<<endl;
    //cerr<<"B = "<<B<<endl;
    //cerr<<"C = "<<C<<endl;
    //cerr<<"alpha = "<<alpha<<endl;
    //cerr<<"Correct result = "<<C2<<endl;
#endif
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.size() == C.rowsize());
    TMVAssert(A.size() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct() == NonConj);
    TMVAssert(!C.isrm());
    TMVAssert(!C.isconj());

    if (SameStorage(A,C) && A.stepj() == C.stepi()) {
      // Then need temporary (see below)
      if (A.isrm()) {
	UpperTriMatrix<Ta,NonUnitDiag,RowMajor> A2 = A;
	ColMultMM<add>(alpha,A2,B,C);
      }
      else {
	UpperTriMatrix<Ta,NonUnitDiag,ColMajor> A2 = A;
	ColMultMM<add>(alpha,A2,B,C);
      }
    } else {
      if (A.isunit()) {
	if (B.isunit()) {
	  T* Cjj = C.ptr();
	  const size_t Cds = C.stepi()+C.stepj();
	  for(size_t j=0,jj=1;j<C.rowsize();++j,++jj,Cjj+=Cds) {
	    // C.col(j) (+=) alpha*A*B.col(j)
	    //
	    // C.col(j) (+=) alpha*A.Cols(j,N)*B.col(j,j,N)
	    //
	    // C(j,j) (+=) alpha*A.row(j,j,N) * B.col(j,j,N)
	    // C.col(j,0,j) (+=) alpha*A.SubMatrix(0,j,j,N)*B.col(j,j,N)
	    // C.col(j,j+1,N) (+=) alpha*A.SubTriMatrix(j,N)*B.col(j,j+1,N)
	    //
	    // Requirements on storage: 
	    //   B can be stored in either triangle
	    //   A cannot be stored in C's lower triangle

	    const size_t N = A.size();

	    T newcjj = A.row(j,jj,N)*B.col(j,jj,N) + T(1);
	    if (alpha != T(1)) newcjj *= alpha;
	    if (add) newcjj += *Cjj;

	    if (add) C.col(j,0,j) += alpha * A.col(j,0,j);
	    else C.col(j,0,j) = alpha * A.col(j,0,j);

	    C.col(j,0,j) += alpha * A.SubMatrix(0,j,jj,N) * B.col(j,jj,N);

	    MultMV<add>(alpha,A.SubTriMatrix(jj,N),B.col(j,jj,N),C.col(j,jj,N));

	    *Cjj = newcjj;
	  }
	} else {
	  const Tb* Bjj = B.cptr();
	  T* Cjj = C.ptr();
	  const size_t Bds = B.stepi()+B.stepj();
	  const size_t Cds = C.stepi()+C.stepj();
	  for(size_t j=0,jj=1;j<C.rowsize();++j,++jj,Bjj+=Bds,Cjj+=Cds) {
	    const size_t N = A.size();

	    T xBjj = B.isconj()?CONJ(*Bjj):*Bjj;
	    T newcjj = A.row(j,jj,N)*B.col(j,jj,N) + xBjj;
	    if (alpha != T(1)) {
	      newcjj *= alpha;
	      xBjj *= alpha; // xBjj is now alpha*B(j,j)
	    }
	    if (add) newcjj += *Cjj;

	    if (add) C.col(j,0,j) += xBjj * A.col(j,0,j);
	    else C.col(j,0,j) = xBjj * A.col(j,0,j);

	    C.col(j,0,j) += alpha * A.SubMatrix(0,j,jj,N) * B.col(j,jj,N);

	    MultMV<add>(alpha,A.SubTriMatrix(jj,N),B.col(j,jj,N),C.col(j,jj,N));

	    *Cjj = newcjj;
	  }
	}
      } else {
	if (B.isunit()) {
	  const Ta* Ajj = A.cptr();
	  T* Cjj = C.ptr();
	  const size_t Ads = A.stepi()+A.stepj();
	  const size_t Cds = C.stepi()+C.stepj();
	  for(size_t j=0,jj=1;j<C.rowsize();++j,++jj,Ajj+=Ads,Cjj+=Cds) {
	    const size_t N = A.size();

	    Ta xAjj = A.isconj() ? CONJ(*Ajj) : *Ajj;
	    T newcjj = A.row(j,jj,N)*B.col(j,jj,N) + xAjj;
	    if (alpha != T(1)) newcjj *= alpha;
	    if (add) newcjj += *Cjj;

	    if (add) C.col(j,0,j) += alpha * A.col(j,0,j);
	    else C.col(j,0,j) = alpha * A.col(j,0,j);

	    C.col(j,0,j) += alpha * A.SubMatrix(0,j,jj,N) * B.col(j,jj,N);

	    MultMV<add>(alpha,A.SubTriMatrix(jj,N),B.col(j,jj,N),C.col(j,jj,N));

	    *Cjj = newcjj;
	  }
	} else {
	  const Tb* Bjj = B.cptr();
	  T* Cjj = C.ptr();
	  const size_t Bds = B.stepi()+B.stepj();
	  const size_t Cds = C.stepi()+C.stepj();
	  for(size_t j=0,jj=1;j<C.rowsize();++j,++jj,Bjj+=Bds,Cjj+=Cds) {
	    const size_t N = A.size();

	    T xBjj = B.isconj() ? CONJ(*Bjj) : *Bjj;
	    T newcjj = A.row(j,j,N)*B.col(j,j,N);
	    if (alpha != T(1)) {
	      newcjj *= alpha;
	      xBjj *= alpha;
	    }
	    if (add) newcjj += *Cjj;
	    
	    if (add) C.col(j,0,j) += xBjj * A.col(j,0,j);
	    else C.col(j,0,j) = xBjj * A.col(j,0,j);

	    C.col(j,0,j) += alpha * A.SubMatrix(0,j,jj,N) * B.col(j,jj,N);

	    MultMV<add>(alpha,A.SubTriMatrix(jj,N),B.col(j,jj,N),C.col(j,jj,N));

	    *Cjj = newcjj;
	  }
	}
      }
    }
#ifdef XDEBUG
    if (Norm(C2-C) > 0.001*(ABS(alpha)*Norm(A0)*Norm(B0)+
	  (add?Norm(C0):RealType(T)(0)))) {
      cerr<<"ColMultMM: alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"Aptr = "<<A.cptr()<<", Bptr = "<<B.cptr();
      cerr<<", Cptr = "<<C.cptr()<<endl;
      cerr<<"--> C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

  template <bool add, class T, class Ta, class Tb> inline void DoMultMM(
      const T alpha, const GenUpperTriMatrix<Ta>& A,
      const GenLowerTriMatrix<Tb>& B, const MatrixView<T>& C)
    // C (+)= alpha * A * B
    // This is designed to work even if A,B are in same storage as C
  {
#ifdef XDEBUG
    //cerr<<"MultMM: "<<alpha<<"  "<<Type(A)<<"  "<<A<<"  "<<Type(B)<<"  "<<B<<endl;
    //cerr<<Type(C)<<"  "<<C<<endl;
    Matrix<T> C2 = C;
    Matrix<T> C0 = C;
    Matrix<Ta> A0 = A;
    Matrix<Tb> B0 = B;
    if (add) C2 += alpha*A0*B0;
    else C2 = alpha*A0*B0;
#endif
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.size() == C.rowsize());
    TMVAssert(A.size() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct() == NonConj);

    const size_t nb = TRI_MM_BLOCKSIZE;
    const size_t N = A.size();

    if (N <= TRI_MM_BLOCKSIZE2) {
      if (C.isrm()) 
	ColMultMM<add>(alpha,B.Transpose(),A.Transpose(),C.Transpose());
      else ColMultMM<add>(alpha,A,B,C);
    } else {
      size_t k = N/2;
      if (k > nb) k = k/nb*nb;

      // [ A00 A01 ] [ B00  0  ] = [ A00 B00 + A01 B10   A01 B11 ]
      // [  0  A11 ] [ B10 B11 ]   [      A11 B10        A11 B11 ]

      ConstUpperTriMatrixView<Ta> A00 = A.SubTriMatrix(0,k);
      ConstMatrixView<Ta> A01 = A.SubMatrix(0,k,k,N);
      ConstUpperTriMatrixView<Ta> A11 = A.SubTriMatrix(k,N);
      ConstLowerTriMatrixView<Tb> B00 = B.SubTriMatrix(0,k);
      ConstMatrixView<Tb> B10 = B.SubMatrix(k,N,0,k);
      ConstLowerTriMatrixView<Tb> B11 = B.SubTriMatrix(k,N);
      MatrixView<T> C00 = C.SubMatrix(0,k,0,k);
      MatrixView<T> C01 = C.SubMatrix(0,k,k,N);
      MatrixView<T> C10 = C.SubMatrix(k,N,0,k);
      MatrixView<T> C11 = C.SubMatrix(k,N,k,N);

      DoMultMM<add>(alpha,A00,B00,C00);
      C00 += alpha*A01*B10;
      if (SameStorage(A01,C10)) {
	if (SameStorage(B10,C01)) {
	  // This is the only case where we need temporary storage, and
	  // I don't image that it is often needed, but it's worth checking.
	  Matrix<T> A01x = A01;
	  MultMM<add>(alpha,A11,B10,C10);
	  MultMM<add>(alpha,B11.Transpose(),A01x.Transpose(),C01.Transpose());
	} else {
	  MultMM<add>(alpha,B11.Transpose(),A01.Transpose(),C01.Transpose());
	  MultMM<add>(alpha,A11,B10,C10);
	}
      } else {
	MultMM<add>(alpha,A11,B10,C10);
	MultMM<add>(alpha,B11.Transpose(),A01.Transpose(),C01.Transpose());
      }
      DoMultMM<add>(alpha,A11,B11,C11);
    }
#ifdef XDEBUG
    if (Norm(C2-C) > 0.001*(ABS(alpha)*Norm(A0)*Norm(B0)+
	  (add?Norm(C0):RealType(T)(0)))) {
      cerr<<"DoMultMM: alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"Aptr = "<<A.cptr()<<", Bptr = "<<B.cptr();
      cerr<<", Cptr = "<<C.cptr()<<endl;
      cerr<<"--> C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

  template <bool add, class T, class Ta, class Tb> void MultMM(
      const T alpha, const GenUpperTriMatrix<Ta>& A,
      const GenLowerTriMatrix<Tb>& B, const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"MultMM: "<<alpha<<"  "<<Type(A)<<"  "<<A<<"  "<<Type(B)<<"  "<<B<<endl;
    //cerr<<Type(C)<<"  "<<C<<endl;
    Matrix<T> C2 = C;
    Matrix<T> C0 = C;
    Matrix<Ta> A0 = A;
    Matrix<Tb> B0 = B;
    if (add) C2 += alpha*A0*B0;
    else C2 = alpha*A0*B0;
#endif
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.size() == C.rowsize());

    const size_t N = A.size();

    if (N==0) return;
    else if (alpha == T(0)) {
      if (!add) C.Zero();
    }
    else if (C.isconj()) 
      DoMultMM<add>(CONJ(alpha),A.Conjugate(),B.Conjugate(),C.Conjugate());
    else
      DoMultMM<add>(alpha,A,B,C);

#ifdef XDEBUG
    if (Norm(C2-C) > 0.001*(ABS(alpha)*Norm(A0)*Norm(B0)+
	  (add?Norm(C0):RealType(T)(0)))) {
      cerr<<"MultMM: alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"Aptr = "<<A.cptr()<<", Bptr = "<<B.cptr();
      cerr<<", Cptr = "<<C.cptr()<<endl;
      cerr<<"--> C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif

  }

  //
  // MultMM: M = L * U
  //

  template <bool add, class T, class Ta, class Tb> inline void ColMultMM(
      const T alpha, const GenLowerTriMatrix<Ta>& A,
      const GenUpperTriMatrix<Tb>& B, const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"MultMM: "<<alpha<<"  "<<Type(A)<<"  "<<A<<"  "<<Type(B)<<"  "<<B<<endl;
    //cerr<<Type(C)<<"  "<<C<<endl;
    Matrix<T> C0 = C;
    Matrix<Ta> A0 = A;
    Matrix<Tb> B0 = B;
    Matrix<T> C2 = C;
    if (add) C2 += alpha*A0*B0;
    else C2 = alpha*A0*B0;
    //cerr<<"Start ColMultMM: L * U\n";
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
    //cerr<<"C = "<<Type(C)<<"  "<<C<<endl;
    //cerr<<"Correct result = "<<C2<<endl;
#endif

    if (SameStorage(A,C) && A.stepi() == C.stepj()) {
      // Then need temporary (see below)
      if (A.isrm()) {
	LowerTriMatrix<Ta,NonUnitDiag,RowMajor> A2 = A;
	ColMultMM<add>(alpha,A2,B,C);
      }
      else {
	LowerTriMatrix<Ta,NonUnitDiag,ColMajor> A2 = A;
	ColMultMM<add>(alpha,A2,B,C);
      }
    } else {
      if (A.isunit()) {
	if (B.isunit()) {
	  const size_t Cds = C.stepi()+C.stepj();
	  T* Cjj = C.ptr()+(C.rowsize()-1)*Cds;
	  for(size_t jj=C.rowsize(),j=jj-1;jj>0;--jj,--j, Cjj-=Cds) {
	    // jj = j+1
	    // C.col(j) (+=) alpha*A*B.col(j)
	    //
	    // C.col(j) (+=) alpha*A.Cols(0,j+1)*B.col(j,0,j+1)
	    //
	    // C(j,j) (+=) alpha*A.row(j,0,j+1)*B.col(j,0,j+1)
	    // C.col(j,j+1,N) (+=) alpha*A.SubMatrix(j+1,N,0,j+1)*B.col(j,0,j+1)
	    // C.col(j,0,j) (+=) alpha*A.SubTriMatrix(0,j)*B.col(j,0,j)
	    //
	    // Requirements on storage: 
	    //   B can be stored in either triangle
	    //   A cannot be stored in C's upper triangle

	    const size_t N = A.size();

	    T newcjj = A.row(j,0,j)*B.col(j,0,j) + T(1);
	    if (alpha != T(1)) newcjj *= alpha;
	    if (add) newcjj += *Cjj;

	    if (add) C.col(j,jj,N) += alpha * A.col(j,jj,N);
	    else C.col(j,jj,N) = alpha * A.col(j,jj,N);

	    C.col(j,jj,N) += alpha * A.SubMatrix(jj,N,0,j) * B.col(j,0,j);

	    MultMV<add>(alpha,A.SubTriMatrix(0,j),B.col(j,0,j),C.col(j,0,j));

	    *Cjj = newcjj;
	  }
	} else {
	  const size_t Nm1 = C.rowsize()-1;
	  const size_t Bds = B.stepi()+B.stepj();
	  const size_t Cds = C.stepi()+C.stepj();
	  const Tb* Bjj = B.cptr()+Nm1*Bds;
	  T* Cjj = C.ptr()+Nm1*Cds;
	  for(size_t jj=C.rowsize(),j=jj-1;jj>0;--jj,--j, Bjj-=Bds,Cjj-=Cds) {

	    const size_t N = A.size();

	    T xBjj = B.isconj()?CONJ(*Bjj):*Bjj;
	    T newcjj = A.row(j,0,j)*B.col(j,0,j) + xBjj;
	    if (alpha != T(1)) {
	      newcjj *= alpha;
	      xBjj *= alpha;
	    }
	    if (add) newcjj += *Cjj;

	    if (add) C.col(j,jj,N) += xBjj * A.col(j,jj,N);
	    else C.col(j,jj,N) = xBjj * A.col(j,jj,N);

	    C.col(j,jj,N) += alpha * A.SubMatrix(jj,N,0,j) * B.col(j,0,j);

	    MultMV<add>(alpha,A.SubTriMatrix(0,j),B.col(j,0,j),C.col(j,0,j));

	    *Cjj = newcjj;
	  }
	}
      } else {
	if (B.isunit()) {
	  const size_t Nm1 = C.rowsize()-1;
	  const size_t Ads = A.stepi()+A.stepj();
	  const size_t Cds = C.stepi()+C.stepj();
	  const Ta* Ajj = A.cptr()+Nm1*Ads;
	  T* Cjj = C.ptr()+Nm1*Cds;
	  for(size_t jj=C.rowsize(),j=jj-1;jj>0;--jj,--j, Ajj-=Ads,Cjj-=Cds) {

	    const size_t N = A.size();

	    T newcjj = A.row(j,0,j)*B.col(j,0,j) + (A.isconj()?CONJ(*Ajj):*Ajj);
	    if (alpha != T(1)) newcjj *= alpha;
	    if (add) newcjj += *Cjj;

	    if (add) C.col(j,jj,N) += alpha * A.col(j,jj,N);
	    else C.col(j,jj,N) = alpha * A.col(j,jj,N);

	    C.col(j,jj,N) += alpha * A.SubMatrix(jj,N,0,j) * B.col(j,0,j);

	    MultMV<add>(alpha,A.SubTriMatrix(0,j),B.col(j,0,j),C.col(j,0,j));

	    *Cjj = newcjj;
	  }
	} else {
	  const size_t Nm1 = C.rowsize()-1;
	  const size_t Bds = B.stepi()+B.stepj();
	  const size_t Cds = C.stepi()+C.stepj();
	  const Tb* Bjj = B.cptr()+Nm1*Bds;
	  T* Cjj = C.ptr()+Nm1*Cds;
	  for(size_t jj=C.rowsize(),j=jj-1;jj>0;--jj,--j, Bjj-=Bds,Cjj-=Cds) {

	    const size_t N = A.size();

	    T xBjj = B.isconj() ? CONJ(*Bjj):*Bjj;
	    T newcjj = A.row(j,0,j+1)*B.col(j,0,j+1);
	    if (alpha != T(1)) {
	      newcjj *= alpha;
	      xBjj *= alpha;
	    }
	    if (add) newcjj += *Cjj;

	    if (add) C.col(j,jj,N) += xBjj * A.col(j,jj,N);
	    else C.col(j,jj,N) = xBjj * A.col(j,jj,N);
	    C.col(j,jj,N) += alpha * A.SubMatrix(jj,N,0,j) * B.col(j,0,j);

	    MultMV<add>(alpha,A.SubTriMatrix(0,j),B.col(j,0,j),C.col(j,0,j));

	    *Cjj = newcjj;
	  }
	}
      }
    }
#ifdef XDEBUG
    if (Norm(C2-C) > 0.001*(ABS(alpha)*Norm(A0)*Norm(B0)+
	  (add?Norm(C0):RealType(T)(0)))) {
      cerr<<"MultMM: alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"Aptr = "<<A.cptr()<<", Bptr = "<<B.cptr();
      cerr<<", Cptr = "<<C.cptr()<<endl;
      cerr<<"--> C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

  template <bool add, class T, class Ta, class Tb> inline void DoMultMM(
      const T alpha, const GenLowerTriMatrix<Ta>& A,
      const GenUpperTriMatrix<Tb>& B, const MatrixView<T>& C)
    // C (+)= alpha * A * B
    // This is designed to work even if A,B are in same storage as C
  {
#ifdef XDEBUG
    //cerr<<"MultMM: "<<alpha<<"  "<<Type(A)<<"  "<<A<<"  "<<Type(B)<<"  "<<B<<endl;
    //cerr<<Type(C)<<"  "<<C<<endl;
    Matrix<T> C0 = C;
    Matrix<Ta> A0 = A;
    Matrix<Tb> B0 = B;
    Matrix<T> C2 = C;
    if (add) C2 += alpha*A0*B0;
    else C2 = alpha*A0*B0;
#endif
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.size() == C.rowsize());
    TMVAssert(A.size() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct() == NonConj);

    const size_t nb = TRI_MM_BLOCKSIZE;
    const size_t N = A.size();

    if (N <= TRI_MM_BLOCKSIZE2) {
      if (C.isrm()) 
	ColMultMM<add>(alpha,B.Transpose(),A.Transpose(),C.Transpose());
      else
	ColMultMM<add>(alpha,A,B,C);
    } else {
      size_t k = N/2;
      if (k > nb) k = k/nb*nb;

      // [ A00  0  ] [ B00 B01 ] = [ A00 B00       A00 B01      ]
      // [ A10 A11 ] [  0  B11 ]   [ A10 B00  A10 B01 + A11 B11 ]

      ConstLowerTriMatrixView<Ta> A00 = A.SubTriMatrix(0,k);
      ConstMatrixView<Ta> A10 = A.SubMatrix(k,N,0,k);
      ConstLowerTriMatrixView<Ta> A11 = A.SubTriMatrix(k,N);
      ConstUpperTriMatrixView<Tb> B00 = B.SubTriMatrix(0,k);
      ConstMatrixView<Tb> B01 = B.SubMatrix(0,k,k,N);
      ConstUpperTriMatrixView<Tb> B11 = B.SubTriMatrix(k,N);
      MatrixView<T> C00 = C.SubMatrix(0,k,0,k);
      MatrixView<T> C01 = C.SubMatrix(0,k,k,N);
      MatrixView<T> C10 = C.SubMatrix(k,N,0,k);
      MatrixView<T> C11 = C.SubMatrix(k,N,k,N);

      DoMultMM<add>(alpha,A11,B11,C11);
      C11 += alpha*A10*B01;
      if (SameStorage(A10,C01)) {
	if (SameStorage(B01,C10)) {
	  // This is the only case where we need temporary storage, and
	  // I don't image that it is often needed, but it's worth checking.
	  Matrix<T> A10x = A10;
	  MultMM<add>(alpha,A00,B01,C01);
	  MultMM<add>(alpha,B00.Transpose(),A10x.Transpose(),C10.Transpose());
	} else {
	  MultMM<add>(alpha,B00.Transpose(),A10.Transpose(),C10.Transpose());
	  MultMM<add>(alpha,A00,B01,C01);
	}
      } else {
	MultMM<add>(alpha,A00,B01,C01);
	MultMM<add>(alpha,B00.Transpose(),A10.Transpose(),C10.Transpose());
      }
      DoMultMM<add>(alpha,A00,B00,C00);
    }
#ifdef XDEBUG
    if (Norm(C2-C) > 0.001*(ABS(alpha)*Norm(A0)*Norm(B0)+
	  (add?Norm(C0):RealType(T)(0)))) {
      cerr<<"MultMM: alpha= "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"Aptr = "<<A.cptr()<<", Bptr = "<<B.cptr();
      cerr<<", Cptr = "<<C.cptr()<<endl;
      cerr<<"--> C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

  template <bool add, class T, class Ta, class Tb> void MultMM(
      const T alpha, const GenLowerTriMatrix<Ta>& A,
      const GenUpperTriMatrix<Tb>& B, const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"MultMM: "<<alpha<<"  "<<Type(A)<<"  "<<A<<"  "<<Type(B)<<"  "<<B<<endl;
    //cerr<<Type(C)<<"  "<<C<<endl;
    Matrix<T> C0 = C;
    Matrix<Ta> A0 = A;
    Matrix<Tb> B0 = B;
    Matrix<T> C2 = C;
    if (add) C2 += alpha*A0*B0;
    else C2 = alpha*A0*B0;
#endif
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.size() == C.rowsize());

    const size_t N = A.size();

    if (N==0) return;
    else if (alpha == T(0)) {
      if (!add) C.Zero();
    }
    else if (C.isconj()) 
      DoMultMM<add>(CONJ(alpha),A.Conjugate(),B.Conjugate(),C.Conjugate());
    else 
      DoMultMM<add>(alpha,A,B,C);

#ifdef XDEBUG
    if (Norm(C2-C) > 0.001*(ABS(alpha)*Norm(A0)*Norm(B0)+
	  (add?Norm(C0):RealType(T)(0)))) {
      cerr<<"MultMM: alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"Aptr = "<<A.cptr()<<", Bptr = "<<B.cptr();
      cerr<<", Cptr = "<<C.cptr()<<endl;
      cerr<<"--> C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_TriMatrixArith_E.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


