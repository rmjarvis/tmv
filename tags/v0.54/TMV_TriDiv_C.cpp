
#include "TMV_Tri.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define TRI_DIV_BLOCKSIZE TMV_BLOCKSIZE
#define TRI_DIV_BLOCKSIZE2 (TMV_BLOCKSIZE/2)
#else
#define TRI_DIV_BLOCKSIZE 64
#define TRI_DIV_BLOCKSIZE2 32
#endif
  
  template <class T, class Ta> inline void RowTriLDivEq(
      const GenUpperTriMatrix<Ta>& A, const UpperTriMatrixView<T>& B)
    // B = A^-1 * B
    // where A is a triangle matrix
  {
    TMVAssert(B.isrm());
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || A.isunit());
    TMVAssert(B.size()>0);
    TMVAssert(B.ct() == NonConj);
    const size_t N = B.size();

    if (A.isunit()) {
      for(int i=N-1; i>=0; --i) 
	B.row(i,i+1,N) -= A.row(i,i+1,N) * B.SubTriMatrix(i+1,N);
    } else {
      TMVAssert(!B.isunit());
      const int Ads = A.stepi()+A.stepj();
      const Ta* Aii = A.cptr() + (N-1)*Ads;
      size_t len=1;
      for(int i=N-1; i>=0; --i,Aii-=Ads,++len) {
	B.row(i,i+1,N) -= A.row(i,i+1,N) * B.SubTriMatrix(i+1,N);
	if (*Aii==Ta(0)) 
	  throw SingularUpperTriMatrix<Ta>(A);
	if (*Aii != Ta(1)) B.row(i,i,N) /= (A.isconj()?CONJ(*Aii):*Aii);
      }
    }
  }

  template <class T, class Ta> inline void ColTriLDivEq(
      const GenUpperTriMatrix<Ta>& A, const UpperTriMatrixView<T>& B)
    // B = A^-1 * B
    // where A is a triangle matrix
  {
    TMVAssert(B.isrm());
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || A.isunit());
    TMVAssert(B.size()>0);
    TMVAssert(B.ct() == NonConj);
    const size_t N = B.size();

    if (A.isunit()) {
      if (B.isunit()) 
	for(int j=N-1; j>=0; --j) {
	  B.SubMatrix(0,j,j+1,N) -= A.col(j,0,j) ^ B.row(j,j+1,N);
	  B.col(j,0,j) -= A.col(j,0,j);
	}
      else
	for(int j=N-1; j>=0; --j) 
	  B.SubMatrix(0,j,j,N) -= A.col(j,0,j) ^ B.row(j,j,N);
    } else {
      TMVAssert(!B.isunit());
      const int Ads = A.stepi()+A.stepj();
      const Ta* Ajj = A.cptr() + (N-1)*Ads;
      size_t len=1;
      for(int j=N-1; j>=0; --j,Ajj-=Ads,++len) {
	if (*Ajj==Ta(0)) 
	  throw SingularUpperTriMatrix<Ta>(A);
	if (*Ajj != Ta(1)) B.row(j,j,N) /= (A.isconj()?CONJ(*Ajj):*Ajj);
	B.SubMatrix(0,j,j,N) -= A.col(j,0,j) ^ B.row(j,j,N);
      }
    } 
  }

  template <class T, class Ta> inline void RowTriLDivEq(
      const GenLowerTriMatrix<Ta>& A, const LowerTriMatrixView<T>& B)
    // B = A^-1 * B
    // where A is a triangle matrix
  {
    TMVAssert(B.isrm());
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || A.isunit());
    TMVAssert(B.size()>0);
    TMVAssert(B.ct() == NonConj);
    const size_t N = B.size();

    if (A.isunit()) {
      for(size_t i=0;i<N;++i) 
	B.row(i,0,i) -= A.row(i,0,i) * B.SubTriMatrix(0,i);
    } else {
      TMVAssert(!B.isunit());
      const int Ads = A.stepi()+A.stepj();
      const Ta* Aii = A.cptr();
      for(size_t i=0;i<N;++i,Aii+=Ads) {
	B.row(i,0,i) -= A.row(i,0,i) * B.SubTriMatrix(0,i);
	if (*Aii==Ta(0)) 
	  throw SingularLowerTriMatrix<Ta>(A);
	if (*Aii != Ta(1)) B.row(i,0,i+1) /= (A.isconj()?CONJ(*Aii):*Aii);
      }
    }
  }

  template <class T, class Ta> inline void ColTriLDivEq(
      const GenLowerTriMatrix<Ta>& A, const LowerTriMatrixView<T>& B)
    // B = A^-1 * B
    // where A is a triangle matrix
  {
    TMVAssert(B.isrm());
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || A.isunit());
    TMVAssert(B.size()>0);
    TMVAssert(B.ct() == NonConj);
    const size_t N = B.size();

    if (A.isunit()) {
      if (B.isunit())
	for(size_t j=0;j<N;++j) {
	  B.col(j,j+1,N) -= A.col(j,j+1,N);
	  B.SubMatrix(j+1,N,0,j) -= A.col(j,j+1,N) ^ B.row(j,0,j);
	}
      else
	for(size_t j=0;j<N;++j) 
	  B.SubMatrix(j+1,N,0,j+1) -= A.col(j,j+1,N) ^ B.row(j,0,j+1);
    } else {
      TMVAssert(!B.isunit());
      const int Ads = A.stepi()+A.stepj();
      const Ta* Ajj = A.cptr();
      for(size_t j=0;j<N;++j,Ajj+=Ads) {
	if (*Ajj==Ta(0)) 
	  throw SingularLowerTriMatrix<Ta>(A);
	if (*Ajj != Ta(1)) B.row(j,0,j+1) /= (A.isconj()?CONJ(*Ajj):*Ajj);
	B.SubMatrix(j+1,N,0,j+1) -= A.col(j,j+1,N) ^ B.row(j,0,j+1);
      }
    }
  }

  template <class T, class Ta> inline void DoTriLDivEq(
      const GenUpperTriMatrix<Ta>& A, const UpperTriMatrixView<T>& B)
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || A.isunit());
    TMVAssert(B.size()>0);
    TMVAssert(B.ct() == NonConj);

    const size_t nb = TRI_DIV_BLOCKSIZE;
    const size_t N = B.size();

    if (N <= TRI_DIV_BLOCKSIZE2) {
      if (B.isrm()) {
	if (A.isrm()) RowTriLDivEq(A,B);
	else ColTriLDivEq(A,B);
      } else {
	if (B.isunit())
	  for(size_t j=0;j<B.rowsize();++j) {
	    B.col(j,0,j) -= A.col(j,0,j);
	    TriLDivEq(A.SubTriMatrix(0,j),B.col(j,0,j));
	  }
	else // B is NonUnitDiag
	  for(size_t j=0;j<B.rowsize();++j)
	    TriLDivEq(A.SubTriMatrix(0,j+1),B.col(j,0,j+1));
      }
    } else {
      size_t k = N/2;
      if (k > nb) k = k/nb*nb;

      ConstUpperTriMatrixView<Ta> A00 = A.SubTriMatrix(0,k);
      ConstMatrixView<Ta> A01 = A.SubMatrix(0,k,k,N);
      ConstUpperTriMatrixView<Ta> A11 = A.SubTriMatrix(k,N);
      UpperTriMatrixView<T> B00 = B.SubTriMatrix(0,k);
      MatrixView<T> B01 = B.SubMatrix(0,k,k,N);
      UpperTriMatrixView<T> B11 = B.SubTriMatrix(k,N);

      DoTriLDivEq(A11,B11);
      B01 -= A01 * B11;
      TriLDivEq(A00,B01);
      DoTriLDivEq(A00,B00);
    }
  }

  template <class T, class Ta> void TriLDivEq(
      const GenUpperTriMatrix<Ta>& A, const UpperTriMatrixView<T>& B)
    // B = A^-1 * B
    // where A is a triangle matrix
  {
#ifdef XDEBUG
    Matrix<T> A0 = A;
    Matrix<T> B0 = B;
#endif
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || A.isunit());

    if (B.isconj()) TriLDivEq(A.Conjugate(),B.Conjugate());
    else DoTriLDivEq(A,B);

#ifdef XDEBUG
    Matrix<T> BB = A0*B;
    if (Norm(BB-B0) > 0.001*Norm(A0)*Norm(B0)) {
      cerr<<"TriLDivEq: Upper/Upper\n";
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"Done: B = "<<B<<endl;
      cerr<<"A*B = "<<BB<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta> inline void DoTriLDivEq(
      const GenLowerTriMatrix<Ta>& A, const LowerTriMatrixView<T>& B)
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || A.isunit());
    TMVAssert(B.size()>0);
    TMVAssert(B.ct() == NonConj);

    const size_t nb = TRI_DIV_BLOCKSIZE;
    const size_t N = B.size();

    if (N <= TRI_DIV_BLOCKSIZE2) {
      if (B.isrm()) {
	if (A.isrm()) RowTriLDivEq(A,B);
	else ColTriLDivEq(A,B);
      } else {
	const size_t N = A.size();
	if (B.isunit())
	  for(size_t j=0;j<B.rowsize();++j) {
	    B.col(j,j+1,N) -= A.col(j,j+1,N);
	    TriLDivEq(A.SubTriMatrix(j+1,N),B.col(j,j+1,N));
	  }
	else // B is NonUnitDiag
	  for(size_t j=0;j<B.rowsize();++j)
	    TriLDivEq(A.SubTriMatrix(j,N),B.col(j,j,N));
      }
    } else {
      size_t k = N/2;
      if (k > nb) k = k/nb*nb;

      ConstLowerTriMatrixView<Ta> A00 = A.SubTriMatrix(0,k);
      ConstMatrixView<Ta> A10 = A.SubMatrix(k,N,0,k);
      ConstLowerTriMatrixView<Ta> A11 = A.SubTriMatrix(k,N);
      LowerTriMatrixView<T> B00 = B.SubTriMatrix(0,k);
      MatrixView<T> B10 = B.SubMatrix(k,N,0,k);
      LowerTriMatrixView<T> B11 = B.SubTriMatrix(k,N);

      DoTriLDivEq(A00,B00);
      B10 -= A10 * B00;
      TriLDivEq(A11,B10);
      TriLDivEq(A11,B11);
    }
  }

  template <class T, class Ta> void TriLDivEq(
      const GenLowerTriMatrix<Ta>& A, const LowerTriMatrixView<T>& B)
    // B = A^-1 * B
    // where A is a triangle matrix
  {
#ifdef XDEBUG
    Matrix<T> A0 = A;
    Matrix<T> B0 = B;
#endif

    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || A.isunit());

    if (B.isconj()) TriLDivEq(A.Conjugate(),B.Conjugate());
    else if (B.size() > 0) DoTriLDivEq(A,B);

#ifdef XDEBUG
    Matrix<T> BB = A0*B;
    if (Norm(BB-B0) > 0.001*Norm(A0)*Norm(B0)) {
      cerr<<"TriLDivEq: Lower/Lower\n";
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"Done: B = "<<B<<endl;
      cerr<<"A*B = "<<BB<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_TriDiv_C.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


