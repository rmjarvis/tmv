
#include "TMV.h"
#include "TMV_Sym.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
  const size_t MM_BLOCKSIZE = TMV_BLOCKSIZE;
#else
  const size_t MM_BLOCKSIZE = 64;
#endif

  //
  // MultMM
  //

  template <class T, class Ta, class Tb> inline void RRowMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"RRowMultMM\n";
#endif
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.uplo() == Lower);

    const size_t N = A.size();
    for(size_t j=0;j<N;++j) {
      // C.row(j) = beta * C.row(j) + alpha * A.row(j,0,j+1) * B.Rows(0,j+1);
      MultMV(alpha,B.Rows(0,j+1).Transpose(),A.row(j,0,j+1),
	    beta,C.row(j));
      C.Rows(0,j) += alpha * A.col(j,0,j) ^ B.row(j);
    }
  }

  template <class T, class Ta, class Tb> inline void CRowMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"CRowMultMM\n";
#endif
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.uplo() == Lower);

    const size_t N = A.size();
    for(int j=N-1;j>=0;--j) {
      // C.row(j) = beta * C.row(j) + alpha * A.row(j,j,N) * B.Rows(j,N);
      MultMV(alpha,B.Rows(j,N).Transpose(),A.row(j,j,N),beta,C.row(j));
      C.Rows(j+1,N) += alpha * A.col(j,j+1,N) ^ B.row(j);
    }
  }

  template <class T, class Ta, class Tb> inline void RowMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    if (A.iscm()) CRowMultMM(alpha,A,B,beta,C);
    else RRowMultMM(alpha,A,B,beta,C);
  }

  template <class T, class Ta, class Tb> inline void ColMultMM(const T alpha,
      const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"ColMultMM\n";
#endif
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.uplo() == Lower);

    for(size_t j=0;j<C.rowsize();++j) 
      // C.col(j) = beta * C.col(j) + alpha * A * B.col(j)
      MultMV(alpha,A,B.col(j),beta,C.col(j));
  }

  template <class T, class Ta, class Tb> inline void NonBlasMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"NonBlockMultMM\n";
#endif
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.uplo() == Lower);

    if (B.isrm() && C.isrm()) RowMultMM(alpha,A,B,beta,C);
    else if (B.iscm() && C.iscm()) ColMultMM(alpha,A,B,beta,C);
    else if (C.colsize() < C.rowsize()) RowMultMM(alpha,A,B,beta,C);
    else ColMultMM(alpha,A,B,beta,C);
  }

#ifdef BLAS
  template <class T, class Ta, class Tb> inline void BlasMultMM(const T alpha,
      const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  { NonBlasMultMM(alpha,A,B,beta,C); }
  template <> inline void BlasMultMM(
      double alpha, const GenSymMatrix<double>& A, const GenMatrix<double>& B,
      double beta, const MatrixView<double>& C)
  {
#ifdef XDEBUG
    //cerr<<"BlasMultMM double\n";
#endif
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != double(0));
    TMVAssert(A.ct()==NonConj);
    TMVAssert(B.ct()==NonConj);
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.isrm() || C.iscm());
    TMVAssert(B.stor() == C.stor());

    if (C.isrm())
      if (A.isrm())
	cblas_dsymm(CblasRowMajor,CblasLeft,CblasLower,
	    C.colsize(),C.rowsize(),alpha,A.cptr(),A.stepi(),
	    B.cptr(),B.stepi(),beta,C.ptr(),C.stepi());
      else
	cblas_dsymm(CblasRowMajor,CblasLeft,CblasUpper,
	    C.colsize(),C.rowsize(),alpha,A.cptr(),A.stepj(),
	    B.cptr(),B.stepi(),beta,C.ptr(),C.stepi());
    else
      if (A.isrm())
	cblas_dsymm(CblasColMajor,CblasLeft,CblasUpper,
	    C.colsize(),C.rowsize(),alpha,A.cptr(),A.stepi(),
	    B.cptr(),B.stepj(),beta,C.ptr(),C.stepj());
      else
	cblas_dsymm(CblasColMajor,CblasLeft,CblasLower,
	    C.colsize(),C.rowsize(),alpha,A.cptr(),A.stepj(),
	    B.cptr(),B.stepj(),beta,C.ptr(),C.stepj());
  }
  template <> inline void BlasMultMM(
      complex<double> alpha, const GenSymMatrix<complex<double> >& A,
      const GenMatrix<complex<double> >& B,
      complex<double> beta, const MatrixView<complex<double> >& C)
  {
#ifdef XDEBUG
    //cerr<<"BlasMultMM c double\n";
#endif
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != double(0));
    TMVAssert(B.ct()==NonConj);
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.isrm() || C.iscm());
    TMVAssert(B.stor() == C.stor());

    if (A.isherm())
      if (C.isrm())
	if (A.isrm())
	  if (A.isconj()) NonBlasMultMM(alpha,A,B,beta,C);
	  else
	    cblas_zhemm(CblasRowMajor,CblasLeft,CblasLower,
		C.colsize(),C.rowsize(),&alpha,A.cptr(),A.stepi(),
		B.cptr(),B.stepi(),&beta,C.ptr(),C.stepi());
	else
	  if (A.isconj())
	    cblas_zhemm(CblasRowMajor,CblasLeft,CblasUpper,
		C.colsize(),C.rowsize(),&alpha,A.cptr(),A.stepj(),
		B.cptr(),B.stepi(),&beta,C.ptr(),C.stepi());
	  else NonBlasMultMM(alpha,A,B,beta,C);
      else
	if (A.isrm())
	  if (A.isconj())
	    cblas_zhemm(CblasColMajor,CblasLeft,CblasUpper,
		C.colsize(),C.rowsize(),&alpha,A.cptr(),A.stepi(),
		B.cptr(),B.stepj(),&beta,C.ptr(),C.stepj());
	  else NonBlasMultMM(alpha,A,B,beta,C);
	else
	  if (A.isconj()) NonBlasMultMM(alpha,A,B,beta,C);
	  else
	    cblas_zhemm(CblasColMajor,CblasLeft,CblasLower,
		C.colsize(),C.rowsize(),&alpha,A.cptr(),A.stepj(),
		B.cptr(),B.stepj(),&beta,C.ptr(),C.stepj());
    else
      if (A.isconj()) NonBlasMultMM(alpha,A,B,beta,C);
      else if (C.isrm())
	if (A.isrm())
	  cblas_zsymm(CblasRowMajor,CblasLeft,CblasLower,
	      C.colsize(),C.rowsize(),&alpha,A.cptr(),A.stepi(),
	      B.cptr(),B.stepi(),&beta,C.ptr(),C.stepi());
	else
	  cblas_zsymm(CblasRowMajor,CblasLeft,CblasUpper,
	      C.colsize(),C.rowsize(),&alpha,A.cptr(),A.stepj(),
	      B.cptr(),B.stepi(),&beta,C.ptr(),C.stepi());
      else
	if (A.isrm())
	  cblas_zsymm(CblasColMajor,CblasLeft,CblasUpper,
	      C.colsize(),C.rowsize(),&alpha,A.cptr(),A.stepi(),
	      B.cptr(),B.stepj(),&beta,C.ptr(),C.stepj());
	else
	  cblas_zsymm(CblasColMajor,CblasLeft,CblasLower,
	      C.colsize(),C.rowsize(),&alpha,A.cptr(),A.stepj(),
	      B.cptr(),B.stepj(),&beta,C.ptr(),C.stepj());
  }
#ifndef NOFLOAT
  template <> inline void BlasMultMM(
      float alpha, const GenSymMatrix<float>& A, const GenMatrix<float>& B,
      float beta, const MatrixView<float>& C)
  {
#ifdef XDEBUG
    //cerr<<"BlasMultMM float\n";
#endif
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != float(0));
    TMVAssert(A.ct()==NonConj);
    TMVAssert(B.ct()==NonConj);
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.isrm() || C.iscm());
    TMVAssert(B.stor() == C.stor());

    if (C.isrm())
      if (A.isrm())
	cblas_ssymm(CblasRowMajor,CblasLeft,CblasLower,
	    C.colsize(),C.rowsize(),alpha,A.cptr(),A.stepi(),
	    B.cptr(),B.stepi(),beta,C.ptr(),C.stepi());
      else
	cblas_ssymm(CblasRowMajor,CblasLeft,CblasUpper,
	    C.colsize(),C.rowsize(),alpha,A.cptr(),A.stepj(),
	    B.cptr(),B.stepi(),beta,C.ptr(),C.stepi());
    else
      if (A.isrm())
	cblas_ssymm(CblasColMajor,CblasLeft,CblasUpper,
	    C.colsize(),C.rowsize(),alpha,A.cptr(),A.stepi(),
	    B.cptr(),B.stepj(),beta,C.ptr(),C.stepj());
      else
	cblas_ssymm(CblasColMajor,CblasLeft,CblasLower,
	    C.colsize(),C.rowsize(),alpha,A.cptr(),A.stepj(),
	    B.cptr(),B.stepj(),beta,C.ptr(),C.stepj());
  }
  template <> inline void BlasMultMM(
      complex<float> alpha, const GenSymMatrix<complex<float> >& A,
      const GenMatrix<complex<float> >& B,
      complex<float> beta, const MatrixView<complex<float> >& C)
  {
#ifdef XDEBUG
    //cerr<<"BlasMultMM c float\n";
#endif
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != float(0));
    TMVAssert(B.ct()==NonConj);
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.isrm() || C.iscm());
    TMVAssert(B.stor() == C.stor());

    if (A.isherm())
      if (C.isrm())
	if (A.isrm())
	  if (A.isconj()) NonBlasMultMM(alpha,A,B,beta,C);
	  else
	    cblas_chemm(CblasRowMajor,CblasLeft,CblasLower,
		C.colsize(),C.rowsize(),&alpha,A.cptr(),A.stepi(),
		B.cptr(),B.stepi(),&beta,C.ptr(),C.stepi());
	else
	  if (A.isconj())
	    cblas_chemm(CblasRowMajor,CblasLeft,CblasUpper,
		C.colsize(),C.rowsize(),&alpha,A.cptr(),A.stepj(),
		B.cptr(),B.stepi(),&beta,C.ptr(),C.stepi());
	  else NonBlasMultMM(alpha,A,B,beta,C);
      else
	if (A.isrm())
	  if (A.isconj())
	    cblas_chemm(CblasColMajor,CblasLeft,CblasUpper,
		C.colsize(),C.rowsize(),&alpha,A.cptr(),A.stepi(),
		B.cptr(),B.stepj(),&beta,C.ptr(),C.stepj());
	  else NonBlasMultMM(alpha,A,B,beta,C);
	else
	  if (A.isconj()) NonBlasMultMM(alpha,A,B,beta,C);
	  else
	    cblas_chemm(CblasColMajor,CblasLeft,CblasLower,
		C.colsize(),C.rowsize(),&alpha,A.cptr(),A.stepj(),
		B.cptr(),B.stepj(),&beta,C.ptr(),C.stepj());
    else
      if (A.isconj()) NonBlasMultMM(alpha,A,B,beta,C);
      else if (C.isrm())
	if (A.isrm())
	  cblas_csymm(CblasRowMajor,CblasLeft,CblasLower,
	      C.colsize(),C.rowsize(),&alpha,A.cptr(),A.stepi(),
	      B.cptr(),B.stepi(),&beta,C.ptr(),C.stepi());
	else
	  cblas_csymm(CblasRowMajor,CblasLeft,CblasUpper,
	      C.colsize(),C.rowsize(),&alpha,A.cptr(),A.stepj(),
	      B.cptr(),B.stepi(),&beta,C.ptr(),C.stepi());
      else
	if (A.isrm())
	  cblas_csymm(CblasColMajor,CblasLeft,CblasUpper,
	      C.colsize(),C.rowsize(),&alpha,A.cptr(),A.stepi(),
	      B.cptr(),B.stepj(),&beta,C.ptr(),C.stepj());
	else
	  cblas_csymm(CblasColMajor,CblasLeft,CblasLower,
	      C.colsize(),C.rowsize(),&alpha,A.cptr(),A.stepj(),
	      B.cptr(),B.stepj(),&beta,C.ptr(),C.stepj());
  }
#endif //NOFLOAT
#endif // BLAS

  template <class T, class Ta, class Tb> void DoMultMM(const T alpha,
      const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));

    if (A.uplo() == Upper) {
      if (A.isherm()) DoMultMM(alpha,A.Adjoint(),B,beta,C);
      else DoMultMM(alpha,A.Transpose(),B,beta,C);
    } else
#ifdef BLAS
      if ( (A.isrm() || A.iscm()) && 
	  (B.isrm() || B.iscm()) && (C.isrm() || C.iscm()) && 
	  (B.ct() == NonConj) && (B.stor() == C.stor()) )
	BlasMultMM(alpha,A,B,beta,C);
      else
#endif
	NonBlasMultMM(alpha,A,B,beta,C);
  }

  template <class T, class Ta, class Tb> void FullTempMultMM(const T alpha,
      const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    if (C.isrm()) {
      Matrix<T,RowMajor> tempC = C;
      DoMultMM(alpha,A,B,beta,tempC.View());
      C = tempC;
    } else {
      Matrix<T,ColMajor> tempC = C;
      DoMultMM(alpha,A,B,beta,tempC.View());
      C = tempC;
    }
  }

  template <class T, class Ta, class Tb> void BlockTempMultMM(const T alpha,
      const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    for(size_t j=0;j<C.rowsize();) {
      size_t j2 = min(C.rowsize(),j+MM_BLOCKSIZE);
      if (C.isrm()) {
	Matrix<T,RowMajor> tempB = B.Cols(j,j2);
	DoMultMM(alpha,A,tempB,beta,C.Cols(j,j2));
      } else {
	Matrix<T,ColMajor> tempB = B.Cols(j,j2);
	DoMultMM(alpha,A,tempB,beta,C.Cols(j,j2));
      }
      j = j2;
    }
  }

  template <class T, class Ta, class Tb> void MultMM(const T alpha,
      const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C
  {
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
#ifdef XDEBUG
    Matrix<T> C2 = alpha*Matrix<T>(A)*B+beta*C;
    Matrix<T> C0 = C;
#endif

    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (C.isconj()) MultMM(CONJ(alpha),A.Conjugate(),B.Conjugate(),
	  CONJ(beta),C.Conjugate());
      else if (alpha == T(0)) 
	C *= beta;
      else if (A.SameStorageAs(C)) 
	FullTempMultMM(alpha,A,B,beta,C);
      else if (B.SameStorageAs(C)) 
	if (C.stepi() == B.stepi() && C.stepj() == B.stepj())
	  BlockTempMultMM(alpha,A,B,beta,C);
	else
	  FullTempMultMM(alpha,A,B,beta,C);
      else DoMultMM(alpha, A, B, beta, C);
    }
#ifdef XDEBUG
    if (Norm(C2-C) > 0.001*Norm(C)) {
      cerr<<"MultMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"--> C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta, class Tb> inline void BlockTempMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenSymMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"NonBlasMultMM\n";
#endif
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == C.rowsize());
    TMVAssert(A.size() > 0);
    TMVAssert(alpha != T(0));

    const size_t N = A.size();

    for(size_t j=0;j<N;) {
      size_t j2 = min(N,j+MM_BLOCKSIZE);
      if (C.isrm()) {
	Matrix<T,RowMajor> tempB(N,j2-j);
	tempB.Rows(0,j) = B.SubMatrix(0,j,j,j2);
	tempB.Rows(j,j2) = B.SubSymMatrix(j,j2);
	tempB.Rows(j2,N) = B.SubMatrix(j2,N,j,j2);
	DoMultMM(alpha,A,tempB.View(),beta,C.Cols(j,j2));
      } else {
	Matrix<T,ColMajor> tempB(N,j2-j);
	tempB.Rows(0,j) = B.SubMatrix(0,j,j,j2);
	tempB.Rows(j,j2) = B.SubSymMatrix(j,j2);
	tempB.Rows(j2,N) = B.SubMatrix(j2,N,j,j2);
	DoMultMM(alpha,A,tempB.View(),beta,C.Cols(j,j2));
      }
      j = j2;
    }
  }

  template <class T, class Ta, class Tb> void FullTempMultMM(const T alpha,
      const GenSymMatrix<Ta>& A, const GenSymMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    if (C.isrm()) {
      Matrix<T,RowMajor> tempC = C;
      BlockTempMultMM(alpha,A,B,beta,tempC.View());
      C = tempC;
    } else {
      Matrix<T,ColMajor> tempC = C;
      BlockTempMultMM(alpha,A,B,beta,tempC.View());
      C = tempC;
    }
  }

  template <class T, class Ta, class Tb> void MultMM(const T alpha,
      const GenSymMatrix<Ta>& A, const GenSymMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == C.rowsize());
#ifdef XDEBUG
    Matrix<T> C2 = alpha*Matrix<T>(A)*Matrix<T>(B)+beta*C;
    Matrix<T> C0 = C;
#endif

    if (A.size() > 0) {
      if (C.isconj()) MultMM(CONJ(alpha),A.Conjugate(),B.Conjugate(),
	  CONJ(beta),C.Conjugate());
      else if (alpha == T(0)) 
	C *= beta;
      else if (A.SameStorageAs(C) || B.SameStorageAs(C))
	FullTempMultMM(alpha,A,B,beta,C);
      else BlockTempMultMM(alpha, A, B, beta, C);
    }

#ifdef XDEBUG
    if (Norm(C2-C) > 0.001*Norm(C)) {
      cerr<<"MultMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"--> C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_SymMatrixArith_D.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


