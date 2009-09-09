
#include "TMV.h"
#include "TMV_Tri.h"

//#define XDEBUG

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

  template <class T, class Ta, class Tb> inline void ColMultMM(
      const T alpha, const GenUpperTriMatrix<Ta>& A,
      const GenLowerTriMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    //cerr<<"ColMultMM UL\n";
    //cerr<<"A = "<<A<<endl;
    //cerr<<"B = "<<B<<endl;
    //cerr<<"C = "<<C<<endl;
    //cerr<<"alpha = "<<alpha<<endl;
    //cerr<<"beta = "<<beta<<endl;
    //cerr<<"Correct result = "<<alpha*Matrix<T>(A)*Matrix<T>(B)+beta*C<<endl;

    if (A.SameStorageAs(C) && A.stepj() == C.stepi()) {
      //cerr<<"need temp\n";
      // Then need temporary (see below)
      if (A.isrm()) {
	UpperTriMatrix<Ta,NonUnitDiag,RowMajor> A2 = A;
	ColMultMM(alpha,A2,B,beta,C);
      }
      else {
	UpperTriMatrix<Ta,NonUnitDiag,ColMajor> A2 = A;
	ColMultMM(alpha,A2,B,beta,C);
      }
    } else {
      //cerr<<"no temp needed\n";
      for(size_t j=0,jj=1;j<C.rowsize();++j,++jj) {
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

	//cerr<<"j = "<<j<<endl;
	const size_t N = A.size();

	T cjj = A.isunit() ?
	  (A.row(j,jj,N)*B.col(j,jj,N) + (B.isunit() ? T(1) : B(j,j))) :
	    B.isunit() ? (A.row(j,jj,N)*B.col(j,jj,N) + A(j,j)) :
	      A.row(j,j,N)*B.col(j,j,N);
	if (alpha != T(1)) cjj *= alpha;
	if (beta != T(0)) {
	  if (beta == T(1)) cjj += C(j,j);
	  else cjj += beta*C(j,j);
	}

	if (B.isunit())
	  AddVV(alpha,A.col(j,0,j),beta,C.col(j,0,j));
	else
	  AddVV(alpha*B(j,j),A.col(j,0,j),beta,C.col(j,0,j));
	MultMV(alpha,A.SubMatrix(0,j,jj,N),B.col(j,jj,N),T(1),C.col(j,0,j));
	//cerr<<"C(0:"<<j<<","<<j<<") => "<<C.col(j,0,j)<<endl;
	
	MultMV(alpha,A.SubTriMatrix(jj,N),B.col(j,jj,N),beta,C.col(j,jj,N));
	//cerr<<"C("<<jj<<":"<<N<<","<<j<<") => "<<C.col(j,jj,N)<<endl;

	C(j,j) = cjj;
	//cerr<<"cjj = "<<cjj<<endl;
      }
    }
  }

  template <class T, class Ta, class Tb> inline void DoMultMM(
      const T alpha, const GenUpperTriMatrix<Ta>& A,
      const GenLowerTriMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C    
    // This is designed to work even if A,B are in same storage as C
  {
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
	ColMultMM(alpha,B.Transpose(),A.Transpose(),beta,C.Transpose());
      else ColMultMM(alpha,A,B,beta,C);
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

      DoMultMM(alpha,A00,B00,beta,C00);
      C00 += alpha*A01*B10;
      if (A01.SameStorageAs(C10)) {
	if (B10.SameStorageAs(C01)) {
	  // This is the only case where we need temporary storage, and
	  // I don't imagine that it is often needed.
	  // But it's worth checking for of course.
	  Matrix<T> A01x = A01;
	  MultMM(alpha,A11,B10,beta,C10);
	  MultMM(alpha,B11.Transpose(),A01x.Transpose(),beta,C01.Transpose());
	} else {
	  MultMM(alpha,B11.Transpose(),A01.Transpose(),beta,C01.Transpose());
	  MultMM(alpha,A11,B10,beta,C10);
	}
      } else {
	MultMM(alpha,A11,B10,beta,C10);
	MultMM(alpha,B11.Transpose(),A01.Transpose(),beta,C01.Transpose());
      }
      DoMultMM(alpha,A11,B11,beta,C11);
    }
  }

  template <class T, class Ta, class Tb> void MultMM(
      const T alpha, const GenUpperTriMatrix<Ta>& A,
      const GenLowerTriMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"MultMM: "<<alpha<<"  "<<Type(A)<<"  "<<A<<"  "<<Type(B)<<"  "<<B<<endl;
    //cerr<<beta<<"  "<<Type(C)<<"  "<<C<<endl;
    Matrix<T> C2 = C;
    Matrix<T> C0 = C;
    Matrix<T> A0 = A;
    Matrix<T> B0 = B;
    C2 *= beta;
    C2 += alpha*A0*B0;
#endif
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.size() == C.rowsize());

    const size_t N = A.size();

    if (N==0) return;
    else if (alpha == T(0)) 
      C *= beta;
    else if (C.isconj()) 
      DoMultMM(CONJ(alpha),A.Conjugate(),B.Conjugate(),CONJ(beta),
	  C.Conjugate());
    else
      DoMultMM(alpha,A,B,beta,C);

#ifdef XDEBUG
    if (Norm(C2-C) > 0.001*(abs(alpha)*Norm(A0)*Norm(B0)+
	  (beta==T(0)?RealType(T)(0):abs(beta)*Norm(C0)))) {
      cerr<<"MultMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"Aptr = "<<A.cptr()<<", Bptr = "<<B.cptr()<<", Cptr = "<<C.cptr()<<endl;
      cerr<<"--> C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif

  }

  template <class T, class Ta, class Tb> inline void ColMultMM(
      const T alpha, const GenLowerTriMatrix<Ta>& A,
      const GenUpperTriMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    //cerr<<"Start ColMultMM: L * U\n";
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
    //cerr<<"C = "<<Type(C)<<"  "<<C<<endl;
    //cerr<<"Correct result = "<<alpha*Matrix<T>(A)*Matrix<T>(B)+beta*C<<endl;

    if (A.SameStorageAs(C) && A.stepi() == C.stepj()) {
      //cerr<<"need temp\n";
      // Then need temporary (see below)
      if (A.isrm()) {
	LowerTriMatrix<Ta,NonUnitDiag,RowMajor> A2 = A;
	ColMultMM(alpha,A2,B,beta,C);
      }
      else {
	LowerTriMatrix<Ta,NonUnitDiag,ColMajor> A2 = A;
	ColMultMM(alpha,A2,B,beta,C);
      }
    } else {
      //cerr<<"no temp needed\n";
      for(size_t jj=C.rowsize(),j=jj-1;jj>0;--jj,--j) { // jj = j+1
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

	T cjj = A.isunit() ?
	  (A.row(j,0,j)*B.col(j,0,j) + (B.isunit() ? T(1) : B(j,j))) :
	    B.isunit() ? (A.row(j,0,j)*B.col(j,0,j) + A(j,j)) :
	      A.row(j,0,j+1)*B.col(j,0,j+1);
	if (alpha != T(1)) cjj *= alpha;
	if (beta != T(0)) {
	  if (beta == T(1)) cjj += C(j,j);
	  else cjj += beta*C(j,j);
	}

	if (B.isunit())
	  AddVV(alpha,A.col(j,jj,N),beta,C.col(j,jj,N));
	else
	  AddVV(alpha*B(j,j),A.col(j,jj,N),beta,C.col(j,jj,N));
	MultMV(alpha,A.SubMatrix(jj,N,0,j),B.col(j,0,j),T(1),C.col(j,jj,N));
	//cerr<<"C("<<jj<<":"<<N<<","<<j<<") => "<<C.col(j,jj,N)<<endl;

	MultMV(alpha,A.SubTriMatrix(0,j),B.col(j,0,j),beta,C.col(j,0,j));
	//cerr<<"C("<<0<<":"<<j<<","<<j<<") => "<<C.col(j,0,j)<<endl;

	C(j,j) = cjj;
	//cerr<<"C("<<j<<","<<j<<") => "<<C(j,j)<<endl;
      }
    }
  }

  template <class T, class Ta, class Tb> inline void DoMultMM(
      const T alpha, const GenLowerTriMatrix<Ta>& A,
      const GenUpperTriMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C    
    // This is designed to work even if A,B are in same storage as C
  {
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
	ColMultMM(alpha,B.Transpose(),A.Transpose(),beta,C.Transpose());
      else
	ColMultMM(alpha,A,B,beta,C);
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

      DoMultMM(alpha,A11,B11,beta,C11);
      C11 += alpha*A10*B01;
      if (A10.SameStorageAs(C01)) {
	if (B01.SameStorageAs(C10)) {
	  // This is the only case where we need temporary storage, and
	  // I don't image that it is often needed.
	  // But it's worth checking for I guess.
	  Matrix<T> A10x = A10;
	  MultMM(alpha,A00,B01,beta,C01);
	  MultMM(alpha,B00.Transpose(),A10x.Transpose(),beta,C10.Transpose());
	} else {
	  MultMM(alpha,B00.Transpose(),A10.Transpose(),beta,C10.Transpose());
	  MultMM(alpha,A00,B01,beta,C01);
	}
      } else {
	MultMM(alpha,A00,B01,beta,C01);
	MultMM(alpha,B00.Transpose(),A10.Transpose(),beta,C10.Transpose());
      }
      DoMultMM(alpha,A00,B00,beta,C00);
    }
  }

  template <class T, class Ta, class Tb> void MultMM(
      const T alpha, const GenLowerTriMatrix<Ta>& A,
      const GenUpperTriMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"MultMM: "<<alpha<<"  "<<Type(A)<<"  "<<A<<"  "<<Type(B)<<"  "<<B<<endl;
    //cerr<<beta<<"  "<<Type(C)<<"  "<<C<<endl;
    Matrix<T> C0 = C;
    Matrix<T> A0 = A;
    Matrix<T> B0 = B;
    Matrix<T> C2 = C;
    C2 *= beta;
    C2 += alpha*A0*B0;
#endif
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.size() == C.rowsize());

    const size_t N = A.size();

    if (N==0) return;
    else if (alpha == T(0)) 
      C *= beta;
    else if (C.isconj()) 
      DoMultMM(CONJ(alpha),A.Conjugate(),B.Conjugate(),CONJ(beta),
	  C.Conjugate());
    else 
      DoMultMM(alpha,A,B,beta,C);

#ifdef XDEBUG
    if (Norm(C2-C) > 0.001*(abs(alpha)*Norm(A0)*Norm(B0)+
	  (beta==T(0)?RealType(T)(0):abs(beta)*Norm(C0)))) {
      cerr<<"MultMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"Aptr = "<<A.cptr()<<", Bptr = "<<B.cptr()<<", Cptr = "<<C.cptr()<<endl;
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


