
#include "TMV.h"
#include "TMV_Tri.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
  const size_t TRI_MM_BLOCKSIZE = TMV_BLOCKSIZE;
  const size_t TRI_MM_BLOCKSIZE2 = TMV_BLOCKSIZE/2;
#else
  const size_t TRI_MM_BLOCKSIZE = 64;
  const size_t TRI_MM_BLOCKSIZE2 = 32;
#endif

  //
  // MultEqMM: M = U * M
  //

  template <class T, class Ta> void RRMultEqMM(T alpha, 
      const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
#ifdef XDEBUG
    //cerr<<"RRMultEq Upper alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.isrm());
    TMVAssert(B.isrm());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(alpha != T(0));
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.colsize()>0);
    // A,B Same storage is ok
    const size_t N = B.colsize();

    if (A.isunit()) {
      for(size_t i=0; i<N; ++i) {
	B.row(i) += A.row(i,i+1,N) * B.Rows(i+1,N);
      }
    }
    else {
      const int Ads = A.stepi()+1;
      const Ta* Aii = A.cptr();
      for(size_t i=0; i<N; ++i,Aii+=Ads) {
	// B.row(i) = (*Aii) * B.row(i) + A.row(i,i+1,N) * B.Rows(i+1,N);
	MultMV(T(1),B.Rows(i+1,N).Transpose(),A.row(i,i+1,N),
	    T(A.isconj()?CONJ(*Aii):*Aii),B.row(i));
      }
    }
    if (alpha != T(1)) B *= alpha;
  }

  template <class T, class Ta> void CRMultEqMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
#ifdef XDEBUG
    //cerr<<"CRMultEq Upper alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.iscm());
    TMVAssert(B.isrm());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(alpha != T(0));
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.colsize()>0);
    // A,B Same storage is ok
    const size_t N = B.colsize();

    if (A.isunit()) {
      for(size_t j=0; j<N; ++j) 
	B.Rows(0,j) += A.col(j,0,j) ^ B.row(j);
    }
    else {
      const int Ads = A.stepi()+A.stepj();
      const Ta* Ajj = A.cptr();
      for(size_t j=0; j<N; ++j,Ajj+=Ads) {
	B.Rows(0,j) += A.col(j,0,j) ^ B.row(j);
	B.row(j) *= (A.isconj()?CONJ(*Ajj):*Ajj);
      }
    }
    B *= alpha;
  }

  template <class T, class Ta> void CMultEqMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
#ifdef XDEBUG
    //cerr<<"CMultEq Upper alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.colsize()>0);
    TMVAssert(alpha != T(0));
    TMVAssert(B.ct() == NonConj);
    TMVAssert(!A.SameStorageAs(B));

    for(size_t j=0;j<B.rowsize();++j) 
      B.col(j) = alpha * A * B.col(j);
  }

  template <class T, class Ta> inline void NonBlasMultEqMM(T alpha, 
      const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
    // B = alpha * A * B
  {
#ifdef XDEBUG
    //cerr<<"NonBlasMultEqMM: "<<alpha<<"  "<<Type(A)<<"  "<<A<<"  "<<Type(B)<<"  "<<B<<endl;
#endif
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(B.ct() == NonConj);

    const size_t nb = TRI_MM_BLOCKSIZE;
    const size_t N = A.size();

    if (N <= TRI_MM_BLOCKSIZE2) {
      if (A.isrm() && B.isrm())
	RRMultEqMM(alpha,A,B);
      else if (A.iscm() && B.isrm())
	CRMultEqMM(alpha,A,B);
      else if (!(A.iscm() || A.isrm()) || A.SameStorageAs(B)) {
	if (A.isunit()) {
	  UpperTriMatrix<T,UnitDiag,ColMajor> AA = A;
	  NonBlasMultEqMM(alpha,AA,B);
	} else {
	  UpperTriMatrix<T,NonUnitDiag,ColMajor> AA = A;
	  NonBlasMultEqMM(alpha,AA,B);
	}
      } 
      else CMultEqMM(alpha,A,B);
    } else {
      size_t k = N/2;
      if (k > nb) k = k/nb*nb;

      // [ A00 A01 ] [ B0 ] = [ A00 B0 + A01 B1 ]
      // [  0  A11 ] [ B1 ]   [      A11 B1     ]

      ConstUpperTriMatrixView<Ta> A00 = A.SubTriMatrix(0,k);
      ConstUpperTriMatrixView<Ta> A11 = A.SubTriMatrix(k,N);
      ConstMatrixView<Ta> A01 = A.SubMatrix(0,k,k,N);
      MatrixView<T> B0 = B.Rows(0,k);
      MatrixView<T> B1 = B.Rows(k,N);

      NonBlasMultEqMM(alpha,A00,B0);
      B0 += alpha * A01 * B1;
      NonBlasMultEqMM(alpha,A11,B1);
    }
  }

#ifdef BLAS
  template <class T, class Ta> inline void BlasMultEqMM(T alpha,
      const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
  { NonBlasMultEqMM(alpha,A,B); }
  inline void BlasMultEqMM(double alpha, const GenUpperTriMatrix<double>& A, 
      const MatrixView<double>& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(alpha != double(0));
    TMVAssert(A.ct() == NonConj);
    TMVAssert(B.ct() == NonConj);
    bool trans = A.stor() != B.stor();
    cblas_dtrmm(B.isrm() ? CblasRowMajor : CblasColMajor, CblasLeft,
	trans ? CblasLower : CblasUpper,
	trans ? CblasTrans : CblasNoTrans,
	A.isunit() ? CblasUnit : CblasNonUnit,
	B.colsize(),B.rowsize(), alpha,
	A.cptr(), A.isrm() ? A.stepi() : A.stepj(),
	B.ptr(), B.isrm() ? B.stepi() : B.stepj()); 
  }
  inline void BlasMultEqMM(complex<double> alpha,
      const GenUpperTriMatrix<complex<double> >& A,
      const MatrixView<complex<double> >& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(alpha != double(0));
    TMVAssert(B.ct() == NonConj);
    bool trans = A.stor() != B.stor();
    if (!trans && A.isconj()) NonBlasMultEqMM(alpha,A,B);
    else {
      CBLAS_ORDER order = B.isrm() ? CblasRowMajor : CblasColMajor;
      CBLAS_UPLO uplo = trans ? CblasLower : CblasUpper;
      CBLAS_TRANSPOSE Atran = trans ?
	A.isconj() ? CblasConjTrans : CblasTrans : CblasNoTrans;
      cblas_ztrmm(order, CblasLeft, uplo, Atran,
	  A.isunit() ? CblasUnit : CblasNonUnit,
	  B.colsize(),B.rowsize(), &alpha,
	  A.cptr(), A.isrm() ? A.stepi() : A.stepj(),
	  B.ptr(), B.isrm() ? B.stepi() : B.stepj()); 
    }
  }
#ifndef NOFLOAT
  inline void BlasMultEqMM(float alpha, const GenUpperTriMatrix<float>& A, 
      const MatrixView<float>& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(alpha != float(0));
    TMVAssert(A.ct() == NonConj);
    TMVAssert(B.ct() == NonConj);
    bool trans = A.stor() != B.stor();
    cblas_strmm(B.isrm() ? CblasRowMajor : CblasColMajor, CblasLeft,
	trans ? CblasLower : CblasUpper,
	trans ? CblasTrans : CblasNoTrans,
	A.isunit() ? CblasUnit : CblasNonUnit,
	B.colsize(),B.rowsize(), alpha,
	A.cptr(), A.isrm() ? A.stepi() : A.stepj(),
	B.ptr(), B.isrm() ? B.stepi() : B.stepj()); 
  }
  inline void BlasMultEqMM(complex<float> alpha,
      const GenUpperTriMatrix<complex<float> >& A,
      const MatrixView<complex<float> >& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(alpha != float(0));
    TMVAssert(B.ct() == NonConj);
    bool trans = A.stor() != B.stor();
    if (!trans && A.isconj()) NonBlasMultEqMM(alpha,A,B);
    else {
      CBLAS_ORDER order = B.isrm() ? CblasRowMajor : CblasColMajor;
      CBLAS_UPLO uplo = trans ? CblasLower : CblasUpper;
      CBLAS_TRANSPOSE Atran = trans ?
	A.isconj() ? CblasConjTrans : CblasTrans : CblasNoTrans;
      cblas_ctrmm(order, CblasLeft, uplo, Atran,
	  A.isunit() ? CblasUnit : CblasNonUnit,
	  B.colsize(),B.rowsize(), &alpha,
	  A.cptr(), A.isrm() ? A.stepi() : A.stepj(),
	  B.ptr(), B.isrm() ? B.stepi() : B.stepj()); 
    }
  }
#endif // NOFLOAT
#endif // BLAS

  template <class T, class Ta> inline void MultEqMM(T alpha, 
      const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
#ifdef XDEBUG
    //cerr<<"MultEqMM: "<<alpha<<"  "<<Type(A)<<"  "<<A<<"  "<<Type(B)<<"  "<<B<<endl;
#endif
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(B.ct() == NonConj);
#ifdef BLAS
    if ( ((A.isrm() && A.stepi()>0) || (A.iscm() && A.stepj()>0)) && 
	((B.isrm() && B.stepi()>0) || (B.iscm() && B.stepj()>0)) &&
	!A.SameStorageAs(B) )
      BlasMultEqMM(alpha,A,B);
    else
#endif
      NonBlasMultEqMM(alpha,A,B);
  }

  //
  // MultEqMM: M = L * M
  //

  template <class T, class Ta> void RRMultEqMM(T alpha, 
      const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
#ifdef XDEBUG
    //cerr<<"RRMultEq Lower alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.isrm());
    TMVAssert(B.isrm());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(alpha != T(0));
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.colsize()>0);
    // A,B Same storage is ok
    const int N = B.colsize();

    if (A.isunit()) {
      for(int i=N-1; i>=0; --i) {
	B.row(i) += A.row(i,0,i) * B.Rows(0,i);
      }
    }
    else {
      const int Ads = A.stepi()+A.stepj();
      const Ta* Aii = A.cptr()+(N-1)*Ads;
      for(int i=N-1; i>=0; --i,Aii-=Ads) {
	// B.row(i) = (*Aii) * B.row(i) + A.row(i,0,i) * B.Rows(0,i);
	MultMV(T(1),B.Rows(0,i).Transpose(),A.row(i,0,i),
	    T(A.isconj()?CONJ(*Aii):*Aii),B.row(i));
      }
    }
    B *= alpha;
  }

  template <class T, class Ta> void CRMultEqMM(
      T alpha, const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
#ifdef XDEBUG
    //cerr<<"CRMultEq Lower alpha = "<<alpha<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
#endif
    TMVAssert(A.iscm());
    TMVAssert(B.isrm());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.colsize()>0);
    TMVAssert(alpha != T(0));
    TMVAssert(B.ct() == NonConj);
    // A,B Same storage is ok
    const int N = B.colsize();

    if (A.isunit()) {
      for(size_t jj=N,j=jj-1; jj>0; --j,--jj) { // jj = j+1
	B.Rows(jj,N) += A.col(j,jj,N) ^ B.row(j);
      }
    }
    else {
      const int Ads = A.stepi()+A.stepj();
      const Ta* Ajj = A.cptr()+(N-1)*Ads;
      for(size_t jj=N,j=jj-1; jj>0; --j,--jj,Ajj-=Ads) { // jj = j+1
	B.Rows(jj,N) += A.col(j,jj,N) ^ B.row(j);
	B.row(j) *= A.isconj()?CONJ(*Ajj):*Ajj;
      }
    }
    B *= alpha;
  }

  template <class T, class Ta> void CMultEqMM(
      T alpha, const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
#ifdef XDEBUG
    //cerr<<"CMultEq Lower alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.colsize()>0);
    TMVAssert(alpha != T(0));
    TMVAssert(B.ct() == NonConj);
    TMVAssert(!A.SameStorageAs(B));

    for(size_t j=0;j<B.rowsize();++j) 
      B.col(j) = alpha * A * B.col(j);
  }

  template <class T, class Ta> inline void NonBlasMultEqMM(T alpha, 
      const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
    // B = alpha * A * B
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(B.ct() == NonConj);

#ifdef XDEBUG
    //cerr<<"NonBlasMultEqMM: "<<alpha<<"  "<<Type(A)<<"  "<<A<<"  "<<Type(B)<<"  "<<B<<endl;
    //Matrix<T> B0x = B;
    //LowerTriMatrix<Ta> A0 = A;
    //Matrix<T> B2 = alpha * Matrix<T>(A) * Matrix<T>(B);
#endif

    const size_t nb = TRI_MM_BLOCKSIZE;
    const size_t N = A.size();

    if (N <= TRI_MM_BLOCKSIZE2) {
      if (A.isrm() && B.isrm())
	RRMultEqMM(alpha,A,B);
      else if (A.iscm() && B.isrm())
	CRMultEqMM(alpha,A,B);
      else if (!(A.iscm() || A.isrm()) || A.SameStorageAs(B)) {
	if (A.isunit()) {
	  LowerTriMatrix<T,UnitDiag,ColMajor> AA = A;
	  NonBlasMultEqMM(alpha,AA,B);
	} else {
	  LowerTriMatrix<T,NonUnitDiag,ColMajor> AA = A;
	  NonBlasMultEqMM(alpha,AA,B);
	}
      } 
      else CMultEqMM(alpha,A,B);
    } else {
      size_t k = N/2;
      if (k > nb) k = k/nb*nb;

      // [ A00  0  ] [ B0 ] = [      A00 B0     ]
      // [ A10 A11 ] [ B1 ]   [ A10 B0 + A11 B1 ]

      ConstLowerTriMatrixView<Ta> A00 = A.SubTriMatrix(0,k);
      ConstLowerTriMatrixView<Ta> A11 = A.SubTriMatrix(k,N);
      ConstMatrixView<Ta> A10 = A.SubMatrix(k,N,0,k);
      MatrixView<T> B0 = B.Rows(0,k);
      MatrixView<T> B1 = B.Rows(k,N);

      NonBlasMultEqMM(alpha,A11,B1);
      B1 += alpha * A10 * B0;
      NonBlasMultEqMM(alpha,A00,B0);
    }
#ifdef XDEBUG
    //if (Norm(B-B2) > 0.001*(Norm(B)+Norm(A))) {
      //cerr<<"MultEqMM alpha = "<<alpha<<endl;
      //cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      //cerr<<"B = "<<Type(B)<<"  "<<B0x<<endl;
      //cerr<<"B = "<<B<<endl;
      //cerr<<"B2 = "<<B2<<endl;
      //cerr<<"Norm(B2-B) = "<<Norm(B2-B)<<endl;
      //abort();
    //}
#endif
  }

#ifdef BLAS
  template <class T, class Ta> inline void BlasMultEqMM(T alpha,
      const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
  { NonBlasMultEqMM(alpha,A,B); }
  inline void BlasMultEqMM(double alpha, const GenLowerTriMatrix<double>& A, 
      const MatrixView<double>& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(alpha != double(0));
    TMVAssert(A.ct() == NonConj);
    TMVAssert(B.ct() == NonConj);
    bool trans = A.stor() != B.stor();
    cblas_dtrmm(B.isrm() ? CblasRowMajor : CblasColMajor, CblasLeft,
	trans ? CblasUpper : CblasLower,
	trans ? CblasTrans : CblasNoTrans,
	A.isunit() ? CblasUnit : CblasNonUnit,
	B.colsize(),B.rowsize(), alpha,
	A.cptr(), A.isrm() ? A.stepi() : A.stepj(),
	B.ptr(), B.isrm() ? B.stepi() : B.stepj()); 
  }
  inline void BlasMultEqMM(complex<double> alpha,
      const GenLowerTriMatrix<complex<double> >& A,
      const MatrixView<complex<double> >& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(alpha != double(0));
    TMVAssert(B.ct() == NonConj);
    bool trans = A.stor() != B.stor();
    if (!trans && A.isconj()) NonBlasMultEqMM(alpha,A,B);
    else {
      CBLAS_ORDER order = B.isrm() ? CblasRowMajor : CblasColMajor;
      CBLAS_UPLO uplo = trans ? CblasUpper : CblasLower;
      CBLAS_TRANSPOSE Atran = trans ?
	A.isconj() ? CblasConjTrans : CblasTrans : CblasNoTrans;
      cblas_ztrmm(order, CblasLeft, uplo, Atran,
	  A.isunit() ? CblasUnit : CblasNonUnit,
	  B.colsize(),B.rowsize(), &alpha,
	  A.cptr(), A.isrm() ? A.stepi() : A.stepj(),
	  B.ptr(), B.isrm() ? B.stepi() : B.stepj()); 
    }
  }
#ifndef NOFLOAT
  inline void BlasMultEqMM(float alpha, const GenLowerTriMatrix<float>& A, 
      const MatrixView<float>& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(alpha != float(0));
    TMVAssert(A.ct() == NonConj);
    TMVAssert(B.ct() == NonConj);
    TMVAssert(B.ct() == NonConj);
    bool trans = A.stor() != B.stor();
    cblas_strmm(B.isrm() ? CblasRowMajor : CblasColMajor, CblasLeft,
	trans ? CblasUpper : CblasLower,
	trans ? CblasTrans : CblasNoTrans,
	A.isunit() ? CblasUnit : CblasNonUnit,
	B.colsize(),B.rowsize(), alpha,
	A.cptr(), A.isrm() ? A.stepi() : A.stepj(),
	B.ptr(), B.isrm() ? B.stepi() : B.stepj()); 
  }
  inline void BlasMultEqMM(complex<float> alpha,
      const GenLowerTriMatrix<complex<float> >& A,
      const MatrixView<complex<float> >& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(alpha != float(0));
    TMVAssert(B.ct() == NonConj);
    bool trans = A.stor() != B.stor();
    if (!trans && A.isconj()) NonBlasMultEqMM(alpha,A,B);
    else {
      CBLAS_ORDER order = B.isrm() ? CblasRowMajor : CblasColMajor;
      CBLAS_UPLO uplo = trans ? CblasUpper : CblasLower;
      CBLAS_TRANSPOSE Atran = trans ?
	A.isconj() ? CblasConjTrans : CblasTrans : CblasNoTrans;
      cblas_ctrmm(order, CblasLeft, uplo, Atran,
	  A.isunit() ? CblasUnit : CblasNonUnit,
	  B.colsize(),B.rowsize(), &alpha,
	  A.cptr(), A.isrm() ? A.stepi() : A.stepj(),
	  B.ptr(), B.isrm() ? B.stepi() : B.stepj()); 
    }
  }
#endif // NOFLOAT
#endif // BLAS

  template <class T, class Ta> inline void MultEqMM(T alpha, 
      const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
#ifdef XDEBUG
    //cerr<<"MultEqMM: "<<alpha<<"  "<<Type(A)<<"  "<<A<<"  "<<Type(B)<<"  "<<B<<endl;
    Matrix<T> B0 = B;
    LowerTriMatrix<Ta> A0 = A;
    Matrix<T> B2 = alpha * Matrix<T>(A) * Matrix<T>(B);
#endif
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(B.ct() == NonConj);
#ifdef BLAS
    if ( !((B.isrm() && B.stepi()>0) || (B.iscm() && B.stepj()>0)) ||
	!A.SameStorageAs(B)) {
      Matrix<T> BB = alpha*B;
      BlasMultEqMM(T(1),A,BB.View());
      B = BB;
    } else if (!((A.isrm() && A.stepi()>0) || (A.iscm() && A.stepj()>0))) {
      if (A.isunit()) {
	LowerTriMatrix<Ta,UnitDiag,RowMajor> AA = A;
	BlasMultEqMM(alpha,AA,B);
      } else if (IMAG(alpha) == RealType(T)(0)) {
	LowerTriMatrix<Ta,UnitDiag,RowMajor> AA = REAL(alpha)*A;
	BlasMultEqMM(T(1),AA,B);
      } else {
	LowerTriMatrix<T,UnitDiag,RowMajor> AA = alpha*A;
	BlasMultEqMM(T(1),AA,B);
      }
    }
    else
      BlasMultEqMM(alpha,A,B);
#else
    NonBlasMultEqMM(alpha,A,B);
#endif

#ifdef XDEBUG
    if (Norm(B-B2) > 0.001*(Norm(B)+Norm(A))) {
      cerr<<"MultEqMM alpha = "<<alpha<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"B = "<<B<<endl;
      cerr<<"B2 = "<<B2<<endl;
      cerr<<"Norm(B2-B) = "<<Norm(B2-B)<<endl;
      abort();
    }
#endif
  }

  //
  // AddMultMM: M += U * M
  //

  template <class T, class Ta, class Tb> void RRAddMultMM(T alpha, 
      const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"RRAddMult Upper alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.isrm());
    TMVAssert(C.isrm());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha != T(0));
    TMVAssert(C.rowsize()>0);
    TMVAssert(C.colsize()>0);
    TMVAssert(!A.isunit());
    TMVAssert(B.isrm() || !C.SameStorageAs(B));
    // A,C Same storage is ok

    const size_t N = C.colsize();

    for(size_t i=0; i<N; ++i) 
      C.row(i) += alpha * A.row(i,i,N) * B.Rows(i,N);
  }

  template <class T, class Ta, class Tb> void CRAddMultMM(T alpha,
      const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"CRAddMult Upper alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.iscm());
    TMVAssert(B.isrm());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha != T(0));
    TMVAssert(C.rowsize()>0);
    TMVAssert(C.colsize()>0);
    TMVAssert(!A.isunit());
    TMVAssert(C.isrm() || !C.SameStorageAs(B));
    TMVAssert(!C.SameStorageAs(A));

    const size_t N = C.colsize();

    if (C.SameStorageAs(B)) {
      for(size_t j=0; j<N; ++j) {
	C.Rows(0,j) += alpha * A.col(j,0,j) ^ B.row(j);
	C.row(j) += alpha * A(j,j) * B.row(j);
      }
    } else {
      for(size_t j=0; j<N; ++j) 
	C.Rows(0,j+1) += alpha * A.col(j,0,j+1) ^ B.row(j);
    }
  }

  template <class T, class Ta, class Tb> void CAddMultMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"CAddMult Upper alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha != T(0));
    TMVAssert(C.rowsize()>0);
    TMVAssert(C.colsize()>0);
    TMVAssert(!A.SameStorageAs(B) || (A.stepi() == B.stepi()))
    TMVAssert(!A.SameStorageAs(C));

    for(size_t j=0;j<C.rowsize();++j) 
      C.col(j) += alpha * A * B.col(j);
  }

  template <class T, class Ta, class Tb> inline void AddMultMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
    // C += alpha * A * B
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha != T(0));
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct() == NonConj);

#ifdef XDEBUG
    UpperTriMatrix<T> A0 = A;
    Matrix<T> B0 = B;
    Matrix<T> C0 = C;
    Matrix<T> C2 = Matrix<T>(C) + alpha*Matrix<T>(A)*Matrix<T>(B);
#endif

    const size_t nb = TRI_MM_BLOCKSIZE;
    const size_t N = A.size();

    if (N <= TRI_MM_BLOCKSIZE2) {
      if (A.isunit() && ((A.isrm()&&C.isrm()) || (A.iscm()&&B.isrm())) ) {
	const size_t N = A.size();
	if (C.SameStorageAs(B) && N > 1) {
	  Matrix<Tb> BB = alpha*B;
	  AddMultMM(T(1),A.OffDiag(),BB.Rows(1,N),C.Rows(0,N-1));
	  C += BB;
	} else {
	  if (N > 1) AddMultMM(alpha,A.OffDiag(),B.Rows(1,N),C.Rows(0,N-1));
	  C += alpha * B;
	}
      }
      else if (A.isrm() && C.isrm())
	RRAddMultMM(alpha,A,B,C);
      else if (A.iscm() && B.isrm())
	CRAddMultMM(alpha,A,B,C);
      else if (!(A.isrm() || A.iscm())) {
	UpperTriMatrix<T,NonUnitDiag,ColMajor> AA = A;
	AddMultMM(alpha,AA,B,C);
      }
      else CAddMultMM(alpha,A,B,C);
    } else {
      size_t k = N/2;
      if (k > nb) k = k/nb*nb;

      ConstUpperTriMatrixView<Ta> A00 = A.SubTriMatrix(0,k);
      ConstUpperTriMatrixView<Ta> A11 = A.SubTriMatrix(k,N);
      ConstMatrixView<Ta> A01 = A.SubMatrix(0,k,k,N);
      ConstMatrixView<Tb> B0 = B.Rows(0,k);
      ConstMatrixView<Tb> B1 = B.Rows(k,N);
      MatrixView<T> C0 = C.Rows(0,k);
      MatrixView<T> C1 = C.Rows(k,N);

      AddMultMM(alpha,A00,B0,C0);
      C0 += alpha * A01 * B1;
      AddMultMM(alpha,A11,B1,C1);
    }
#ifdef XDEBUG
    if (Norm(C-C2) > 0.001*(Norm(C)+Norm(A)+Norm(B))) {
      cerr<<"AddMultMM alpha = "<<alpha<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"Aptr = "<<A.cptr()<<", Bptr = "<<B.cptr()<<", Cptr = "<<C.cptr()<<endl;
      cerr<<"C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      cerr<<"Norm(C2-C) = "<<Norm(C2-C)<<endl;
      abort();
    }
#endif
  }

  //
  // AddMultMM: M += L * M
  //

  template <class T, class Ta, class Tb> void RRAddMultMM(T alpha, 
      const GenLowerTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"RRAddMult Lower alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.isrm());
    TMVAssert(C.isrm());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha != T(0));
    TMVAssert(C.rowsize()>0);
    TMVAssert(C.colsize()>0);

    const size_t N = C.colsize();

    for(size_t i=0; i<N; ++i) 
      C.row(i) += alpha * A.row(i,0,i+1) * B.Rows(0,i+1);
  }

  template <class T, class Ta, class Tb> void CRAddMultMM(
      T alpha, const GenLowerTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"CRAddMult Lower alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.iscm());
    TMVAssert(B.isrm());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha != T(0));
    TMVAssert(C.rowsize()>0);
    TMVAssert(C.colsize()>0);
    TMVAssert(!A.isunit());

    const size_t N = C.colsize();

    for(size_t j=0; j<N; ++j) 
      C.Rows(j,N) += alpha * A.col(j,j,N) ^ B.row(j);
  }

  template <class T, class Ta, class Tb> void CAddMultMM(
      T alpha, const GenLowerTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"CAddMult Lower alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha != T(0));
    TMVAssert(C.rowsize()>0);
    TMVAssert(C.colsize()>0);
    TMVAssert(!A.SameStorageAs(B) || (A.stepi() == B.stepi()))
    TMVAssert(!A.SameStorageAs(C));

    for(size_t j=0;j<C.rowsize();++j) 
      C.col(j) += alpha * A * B.col(j);
  }

  template <class T, class Ta, class Tb> inline void AddMultMM(
      T alpha, const GenLowerTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
    // C += alpha * A * B
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha != T(0));
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct() == NonConj);

    const size_t nb = TRI_MM_BLOCKSIZE;
    const size_t N = A.size();

    if (N <= TRI_MM_BLOCKSIZE2) {
      if (A.isunit() && ((A.isrm()&&C.isrm()) || (A.iscm()&&B.isrm())) ) {
	const size_t N = A.size();
	if (C.SameStorageAs(B) && N > 1) {
	  Matrix<Tb> BB = alpha*B;
	  AddMultMM(T(1),A.OffDiag(),BB.Rows(0,N-1),C.Rows(1,N));
	  C += BB;
	} else {
	  if (N > 1) AddMultMM(alpha,A.OffDiag(),B.Rows(0,N-1),C.Rows(1,N));
	  C += alpha * B;
	}
      }
      else if (A.isrm() && C.isrm())
	RRAddMultMM(alpha,A,B,C);
      else if (A.iscm() && B.isrm())
	CRAddMultMM(alpha,A,B,C);
      else if (!(A.isrm() || A.iscm())) {
	LowerTriMatrix<T,NonUnitDiag,ColMajor> AA = A;
	AddMultMM(alpha,AA,B,C);
      }
      else CAddMultMM(alpha,A,B,C);
    } else {
      size_t k = N/2;
      if (k > nb) k = k/nb*nb;

      ConstLowerTriMatrixView<Ta> A00 = A.SubTriMatrix(0,k);
      ConstMatrixView<Ta> A10 = A.SubMatrix(k,N,0,k);
      ConstLowerTriMatrixView<Ta> A11 = A.SubTriMatrix(k,N);
      ConstMatrixView<Tb> B0 = B.Rows(0,k);
      ConstMatrixView<Tb> B1 = B.Rows(k,N);
      MatrixView<T> C0 = C.Rows(0,k);
      MatrixView<T> C1 = C.Rows(k,N);

      AddMultMM(alpha,A11,B1,C1);
      C1 += alpha * A10 * B0;
      AddMultMM(alpha,A00,B0,C0);
    }
  }

  template <class T, class Ta, class Tb> void BlockTempMultMM(const T alpha,
      const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C    
  { 
    for (size_t j=0;j<C.rowsize();) {
      size_t j2 = min(C.rowsize(),j+TRI_MM_BLOCKSIZE);
      if (B.isrm()) {
	Matrix<T,RowMajor> B2 = alpha * B.Cols(j,j2);
	MultEqMM(T(1),A,B2.View());
	C.Cols(j,j2) = B2 + beta*C.Cols(j,j2);
      } else {
	Matrix<T,ColMajor> B2 = alpha * B.Cols(j,j2);
	MultEqMM(T(1),A,B2.View());
	C.Cols(j,j2) = B2 + beta*C.Cols(j,j2);
      }
      j = j2;
    }
  }

  template <class T, class Ta, class Tb> void FullTempMultMM(const T alpha,
      const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C    
  { 
    if (B.isrm()) {
      Matrix<T,RowMajor> B2 = alpha * B;
      MultEqMM(T(1),A,B2.View());
      C = B2 + beta*C;
    } else {
      Matrix<T,ColMajor> B2 = alpha * B;
      MultEqMM(T(1),A,B2.View());
      C = B2 + beta*C;
    }
  }

  template <class T, class Ta, class Tb> void MultMM(const T alpha,
      const GenUpperTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C    
  { 
#ifdef XDEBUG
    //cerr<<"MultMM: "<<alpha<<"  "<<Type(A)<<"  "<<A<<"  "<<Type(B)<<"  "<<B<<endl;
    //cerr<<beta<<"  "<<Type(C)<<"  "<<C<<endl;
    Matrix<T> C2 = beta*Matrix<T>(C) + alpha*Matrix<T>(A)*Matrix<T>(B);
    Matrix<T> C0 = C;
#endif
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());

    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (C.isconj()) MultMM(CONJ(alpha),A.Conjugate(),B.Conjugate(),
	  CONJ(beta),C.Conjugate());
      else if (alpha==T(0)) 
	C *= beta;
      else if (C.rowsize() == 1)
	MultMV(alpha,A,B.col(0),beta,C.col(0));
      else if (A.SameStorageAs(C)) 
	FullTempMultMM(alpha,A,B,beta,C);
      else if (beta == T(0)) {
	C = alpha * B;
	MultEqMM(T(1),A,C);
      } 
      else if (C.SameStorageAs(B)) 
	if (C.stepi() == B.stepi() && C.stepj() == B.stepj())
	  BlockTempMultMM(alpha,A,B,beta,C);
	else
	  FullTempMultMM(alpha,A,B,beta,C);
      else {
	C *= beta;
	AddMultMM(alpha,A,B,C);
      } 
    }
#ifdef XDEBUG
    if (Norm(C-C2) > 0.001*(Norm(C)+Norm(A)+Norm(B))) {
      cerr<<"MultMM alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      cerr<<"Norm(C2-C) = "<<Norm(C2-C)<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta, class Tb> void FullTempMultMM(const T alpha,
      const GenLowerTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C    
  { 
    if (B.isrm()) {
      Matrix<T,RowMajor> B2 = alpha * B;
      MultEqMM(T(1),A,B2.View());
      C = B2 + beta*C;
    } else {
      Matrix<T,ColMajor> B2 = alpha * B;
      MultEqMM(T(1),A,B2.View());
      C = B2 + beta*C;
    }
  }

  template <class T, class Ta, class Tb> void BlockTempMultMM(const T alpha,
      const GenLowerTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C    
  { 
    for (size_t j=0;j<C.rowsize();) {
      size_t j2 = min(C.rowsize(),j+TRI_MM_BLOCKSIZE);
      if (B.isrm()) {
	Matrix<T,RowMajor> B2 = alpha * B.Cols(j,j2);
	MultEqMM(T(1),A,B2.View());
	C.Cols(j,j2) = B2 + beta * C.Cols(j,j2);
      } else {
	Matrix<T,ColMajor> B2 = alpha * B.Cols(j,j2);
	MultEqMM(T(1),A,B2.View());
	C.Cols(j,j2) = B2 + beta * C.Cols(j,j2);
      }
      j = j2;
    }
  }

  template <class T, class Ta, class Tb> void MultMM(const T alpha,
      const GenLowerTriMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C    
  { 
#ifdef XDEBUG
    //cerr<<"MultMM: "<<alpha<<"  "<<Type(A)<<"  "<<A<<"  "<<Type(B)<<"  "<<B<<endl;
    //cerr<<beta<<"  "<<Type(C)<<"  "<<C<<endl;
    Matrix<T> C0 = C;
    Matrix<T> A0 = A;
    Matrix<T> B0 = B;
    Matrix<T> C2 = beta*Matrix<T>(C) + alpha*Matrix<T>(A)*Matrix<T>(B);
#endif
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());

    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (C.isconj()) MultMM(CONJ(alpha),A.Conjugate(),B.Conjugate(),
	  CONJ(beta),C.Conjugate());
      else if (alpha==T(0)) 
	C *= beta;
      else if (C.rowsize() == 1)
	MultMV(alpha,A,B.col(0),beta,C.col(0));
      else if (A.SameStorageAs(C)) 
	FullTempMultMM(alpha,A,B,beta,C);
      else if (beta == T(0)) {
	C = alpha * B;
	MultEqMM(T(1),A,C);
      } 
      else if (C.SameStorageAs(B)) 
	if (C.stepi() == B.stepi() && C.stepj() == B.stepj())
	  BlockTempMultMM(alpha,A,B,beta,C);
	else
	  FullTempMultMM(alpha,A,B,beta,C);
      else {
	C *= beta;
	AddMultMM(alpha,A,B,C);
      }
    }
#ifdef XDEBUG
    if (Norm(C-C2) > 0.001*(Norm(C)+Norm(A)+Norm(B))) {
      cerr<<"MultMM alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"Aptr = "<<A.cptr()<<", Bptr = "<<B.cptr();
      cerr<<", Cptr = "<<C.cptr()<<endl;
      cerr<<"C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      cerr<<"Norm(C2-C) = "<<Norm(C2-C)<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_TriMatrixArith_D.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


