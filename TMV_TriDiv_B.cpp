
#include "TMV_Tri.h"
#include "TMV_Diag.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
  const size_t TRI_DIV_BLOCKSIZE = TMV_BLOCKSIZE;
  const size_t TRI_DIV_BLOCKSIZE2 = TMV_BLOCKSIZE/2;
#else
  const size_t TRI_DIV_BLOCKSIZE = 64;
  const size_t TRI_DIV_BLOCKSIZE2 = 32;
#endif
  
  //
  // TriLDivEq M
  //

  template <class T, class Ta> void RowTriLDivEq(
      const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
    // Solve A X = B  where A is an upper triangle matrix
    const size_t N = A.size();
    TMVAssert(B.isrm());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.colsize()>0);
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.ct() == NonConj);

    if (A.isunit()) {
      for(int i=N-1; i>=0; --i) 
	B.row(i) -= A.row(i,i+1,N) * B.Rows(i+1,N);
    } else {
      const int Ads = A.stepi() + A.stepj();
      const Ta* Aii = A.cptr() + (N-1) * Ads;
      for(int i=N-1; i>=0; --i,Aii-=Ads) {
	B.row(i) -= A.row(i,i+1,N) * B.Rows(i+1,N);
	if (*Aii==Ta(0)) tmv_error("Singular Matrix found in UpperTriLDivEq");
	B.row(i) /= (A.isconj() ? CONJ(*Aii) : *Aii);
      }
    }
  }

  template <class T, class Ta> void ColTriLDivEq(
      const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
    // Solve A X = B  where A is an upper triangle matrix
    const size_t N = A.size();
    TMVAssert(B.isrm());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.colsize()>0);
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.ct() == NonConj);

    if (A.isunit()) {
      for(int j=N-1; j>0; --j) 
	B.Rows(0,j) -= A.col(j,0,j) ^ B.row(j);
    } else {
      const int Ads = A.stepi()+A.stepj();
      const Ta* Ajj = A.cptr() + (N-1)*Ads;
      for(int j=N-1; j>=0; --j,Ajj-=Ads) {
	if (*Ajj==Ta(0)) tmv_error("Singular Matrix found in UpperTriLDivEq");
	B.row(j) /= (A.isconj() ? CONJ(*Ajj) : *Ajj);
	B.Rows(0,j) -= A.col(j,0,j) ^ B.row(j);
      }
    } 
  }

  template <class T, class Ta> void RowTriLDivEq(
      const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
    // Solve A X = B  where A is a lower triangle matrix
    const size_t N = A.size();
    TMVAssert(B.isrm());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.colsize()>0);
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.ct() == NonConj);

    if (A.isunit()) {
      for(size_t i=0;i<N;++i) 
	B.row(i) -= A.row(i,0,i) * B.Rows(0,i);
    } else {
      const int Ads = A.stepi()+A.stepj();
      const Ta* Aii = A.cptr();
      for(size_t i=0;i<N;++i,Aii+=Ads) {
	B.row(i) -= A.row(i,0,i) * B.Rows(0,i);
	if (*Aii==Ta(0)) tmv_error("Singular Matrix found in LowerTriLDivEq");
	B.row(i) /= (A.isconj() ? CONJ(*Aii) : *Aii);
      }
    }
  }

  template <class T, class Ta> void ColTriLDivEq(
      const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
    // Solve A X = B  where A is a lower triangle matrix
    const size_t N = A.size();
    TMVAssert(B.isrm());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.colsize()>0);
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.ct() == NonConj);

    if (A.isunit()) {
      for(size_t j=0;j<N;++j) 
	B.Rows(j+1,N) -= A.col(j,j+1,N) ^ B.row(j);
    } else {
      const int Ads = A.stepi()+A.stepj();
      const Ta* Ajj = A.cptr();
      for(size_t j=0;j<N;++j,Ajj+=Ads) {
	if (*Ajj==Ta(0)) tmv_error("Singular Matrix found in LowerTriLDivEq");
	B.row(j) /= (A.isconj() ? CONJ(*Ajj) : *Ajj);
	B.Rows(j+1,N) -= A.col(j,j+1,N) ^ B.row(j);
      }
    } 
  }

  template <class T, class Ta> inline void NonBlasTriLDivEq(
      const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
    // B = A^-1 * B
    // where A is a triangle matrix
  {
    //cerr<<"Upper Tri LDivEq Matrix\n";
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.colsize()>0);
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.ct() == NonConj);

    const size_t nb = TRI_DIV_BLOCKSIZE;
    const size_t N = A.size();

    if (N <= TRI_DIV_BLOCKSIZE2) {
      if (B.isrm()) {
	if (A.isrm()) RowTriLDivEq(A,B);
	else ColTriLDivEq(A,B);
      } else {
	for(size_t j=0;j<B.rowsize();++j) TriLDivEq(A,B.col(j));
      }
    } else {
      size_t k = N/2;
      if (k > nb) k = k/nb*nb;

      ConstUpperTriMatrixView<Ta> A00 = A.SubTriMatrix(0,k);
      ConstUpperTriMatrixView<Ta> A11 = A.SubTriMatrix(k,N);
      ConstMatrixView<Ta> A01 = A.SubMatrix(0,k,k,N);
      MatrixView<T> B0 = B.Rows(0,k);
      MatrixView<T> B1 = B.Rows(k,N);

      NonBlasTriLDivEq(A11,B1);
      B0 -= A01*B1;
      NonBlasTriLDivEq(A00,B0);
    }
  }

  template <class T, class Ta> inline void NonBlasTriLDivEq(
      const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
    // B = A^-1 * B
    // where A is a triangle matrix
  {
    //cerr<<"Lower Tri LDivEq Matrix\n";
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.colsize()>0);
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.ct() == NonConj);

    const size_t nb = TRI_DIV_BLOCKSIZE;
    const size_t N = A.size();

    if (N <= TRI_DIV_BLOCKSIZE2) {
      if (B.isrm()) {
	if (A.isrm()) RowTriLDivEq(A,B);
	else ColTriLDivEq(A,B);
      } else {
	for(size_t j=0;j<B.rowsize();++j) TriLDivEq(A,B.col(j));
      }
    } else {
      size_t k = N/2;
      if (k > nb) k = k/nb*nb;

      ConstLowerTriMatrixView<Ta> A00 = A.SubTriMatrix(0,k);
      ConstLowerTriMatrixView<Ta> A11 = A.SubTriMatrix(k,N);
      ConstMatrixView<Ta> A10 = A.SubMatrix(k,N,0,k);
      MatrixView<T> B0 = B.Rows(0,k);
      MatrixView<T> B1 = B.Rows(k,N);

      NonBlasTriLDivEq(A00,B0);
      B1 -= A10*B0;
      NonBlasTriLDivEq(A11,B1);
    }
  }

#ifdef BLAS
  template <class T, class Ta> inline void BlasTriLDivEq(
      const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
  { NonBlasTriLDivEq(A,B); }
  template <class T, class Ta> inline void BlasTriLDivEq(
      const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
  { NonBlasTriLDivEq(A,B); }
  inline void BlasTriLDivEq(
      const GenUpperTriMatrix<double>& A, const MatrixView<double>& B)
  { 
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.colsize()>0);
    TMVAssert(B.rowsize()>0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(B.ct() == NonConj);
    bool trans = (A.stor() != B.stor());
    double alpha = 1.;
    cblas_dtrsm(B.isrm() ? CblasRowMajor : CblasColMajor,
	CblasLeft, trans ? CblasLower : CblasUpper,
	trans ? CblasTrans : CblasNoTrans,
	A.isunit() ? CblasUnit : CblasNonUnit,
	B.colsize(),B.rowsize(),alpha,
	A.cptr(),A.isrm()?A.stepi():A.stepj(),
	B.ptr(),B.isrm()?B.stepi():B.stepj());
  }
  inline void BlasTriLDivEq(
      const GenLowerTriMatrix<double>& A, const MatrixView<double>& B)
  { 
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.colsize()>0);
    TMVAssert(B.rowsize()>0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(B.ct() == NonConj);
    bool trans = (A.stor() != B.stor());
    double alpha = 1.;
    cblas_dtrsm(B.isrm() ? CblasRowMajor : CblasColMajor,
	CblasLeft, trans ? CblasUpper : CblasLower,
	trans ? CblasTrans : CblasNoTrans,
	A.isunit() ? CblasUnit : CblasNonUnit,
	B.colsize(),B.rowsize(),alpha,
	A.cptr(),A.isrm()?A.stepi():A.stepj(),
	B.ptr(),B.isrm()?B.stepi():B.stepj());
  }
  inline void BlasTriLDivEq(const GenUpperTriMatrix<complex<double> >& A,
      const MatrixView<complex<double> >& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.colsize()>0);
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.ct() == NonConj);
    bool trans = (A.stor() != B.stor());
    if (!trans && A.isconj()) NonBlasTriLDivEq(A,B);
    else {
      CBLAS_ORDER order = B.isrm() ? CblasRowMajor : CblasColMajor;
      CBLAS_UPLO uplo = trans ? CblasLower : CblasUpper;
      CBLAS_TRANSPOSE Atran = trans ? 
	A.isconj() ? CblasConjTrans : CblasTrans : CblasNoTrans;
      complex<double> calpha = 1.;
      cblas_ztrsm(order, CblasLeft, uplo, Atran,
	  A.isunit() ? CblasUnit : CblasNonUnit,
	  B.colsize(),B.rowsize(),&calpha,
	  A.cptr(),A.isrm()?A.stepi():A.stepj(),
	  B.ptr(),B.isrm()?B.stepi():B.stepj());
    }
  }
  inline void BlasTriLDivEq(const GenLowerTriMatrix<complex<double> >& A,
      const MatrixView<complex<double> >& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.colsize()>0);
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.ct() == NonConj);
    bool trans = (A.stor() != B.stor());
    if (!trans && A.isconj()) NonBlasTriLDivEq(A,B);
    else {
      CBLAS_ORDER order = B.isrm() ? CblasRowMajor : CblasColMajor;
      CBLAS_UPLO uplo = trans ? CblasUpper : CblasLower;
      CBLAS_TRANSPOSE Atran = trans ? 
	A.isconj() ? CblasConjTrans : CblasTrans : CblasNoTrans;
      complex<double> calpha = 1.;
      cblas_ztrsm(order, CblasLeft, uplo, Atran,
	  A.isunit() ? CblasUnit : CblasNonUnit,
	  B.colsize(),B.rowsize(),&calpha,
	  A.cptr(),A.isrm()?A.stepi():A.stepj(),
	  B.ptr(),B.isrm()?B.stepi():B.stepj());
    }
  }
#ifndef NOFLOAT
#ifndef MKL
  inline void BlasTriLDivEq(
      const GenUpperTriMatrix<float>& A, const MatrixView<float>& B)
  { 
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.colsize()>0);
    TMVAssert(B.rowsize()>0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(B.ct() == NonConj);
    bool trans = (A.stor() != B.stor());
    float alpha = 1.;
    cblas_strsm(B.isrm() ? CblasRowMajor : CblasColMajor,
	CblasLeft, trans ? CblasLower : CblasUpper,
	trans ? CblasTrans : CblasNoTrans,
	A.isunit() ? CblasUnit : CblasNonUnit,
	B.colsize(),B.rowsize(),alpha,
	A.cptr(),A.isrm()?A.stepi():A.stepj(),
	B.ptr(),B.isrm()?B.stepi():B.stepj());
  }
  inline void BlasTriLDivEq(
      const GenLowerTriMatrix<float>& A, const MatrixView<float>& B)
  { 
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.colsize()>0);
    TMVAssert(B.rowsize()>0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(B.ct() == NonConj);
    bool trans = (A.stor() != B.stor());
    float alpha = 1.;
    cblas_strsm(B.isrm() ? CblasRowMajor : CblasColMajor,
	CblasLeft, trans ? CblasUpper : CblasLower,
	trans ? CblasTrans : CblasNoTrans,
	A.isunit() ? CblasUnit : CblasNonUnit,
	B.colsize(),B.rowsize(),alpha,
	A.cptr(),A.isrm()?A.stepi():A.stepj(),
	B.ptr(),B.isrm()?B.stepi():B.stepj());
  }
  inline void BlasTriLDivEq(
      const GenUpperTriMatrix<complex<float> >& A,
      const MatrixView<complex<float> >& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.colsize()>0);
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.ct() == NonConj);
    bool trans = (A.stor() != B.stor());
    if (!trans && A.isconj()) NonBlasTriLDivEq(A,B);
    else {
      CBLAS_ORDER order = B.isrm() ? CblasRowMajor : CblasColMajor;
      CBLAS_UPLO uplo = trans ? CblasLower : CblasUpper;
      CBLAS_TRANSPOSE Atran = trans ? 
	A.isconj() ? CblasConjTrans : CblasTrans : CblasNoTrans;
      complex<float> calpha = 1.;
      cblas_ctrsm(order, CblasLeft, uplo, Atran,
	  A.isunit() ? CblasUnit : CblasNonUnit,
	  B.colsize(),B.rowsize(),&calpha,
	  A.cptr(),A.isrm()?A.stepi():A.stepj(),
	  B.ptr(),B.isrm()?B.stepi():B.stepj());
    }
  }
  inline void BlasTriLDivEq(
      const GenLowerTriMatrix<complex<float> >& A,
      const MatrixView<complex<float> >& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.colsize()>0);
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.ct() == NonConj);
    bool trans = (A.stor() != B.stor());
    if (!trans && A.isconj()) NonBlasTriLDivEq(A,B);
    else {
      CBLAS_ORDER order = B.isrm() ? CblasRowMajor : CblasColMajor;
      CBLAS_UPLO uplo = trans ? CblasUpper : CblasLower;
      CBLAS_TRANSPOSE Atran = trans ? 
	A.isconj() ? CblasConjTrans : CblasTrans : CblasNoTrans;
      complex<float> calpha = 1.;
      cblas_ctrsm(order, CblasLeft, uplo, Atran,
	  A.isunit() ? CblasUnit : CblasNonUnit,
	  B.colsize(),B.rowsize(),&calpha,
	  A.cptr(),A.isrm()?A.stepi():A.stepj(),
	  B.ptr(),B.isrm()?B.stepi():B.stepj());
    }
  }
#endif // MKL
#endif // FLOAT
#endif // BLAS
  template <class T, class Ta> inline void TriLDivEq(
      const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
#ifdef XDEBUG
    Matrix<T> B0 = B;
#endif

    TMVAssert(A.size() == B.colsize());
    if (B.colsize() > 0 && B.rowsize() > 0) {
      if (B.isconj()) TriLDivEq(A.Conjugate(),B.Conjugate());
      else if (B.rowsize() == 1) TriLDivEq(A,B.col(0));
      else if (B.SameStorageAs(A)) {
	if (A.dt() == NonUnitDiag) {
	  if (A.isrm()) {
	    UpperTriMatrix<Ta,NonUnitDiag,RowMajor> tempA = A;
	    TriLDivEq(tempA,B);
	  } else {
	    UpperTriMatrix<Ta,NonUnitDiag,ColMajor> tempA = A;
	    TriLDivEq(tempA,B);
	  }
	} else {
	  if (A.isrm()) {
	    UpperTriMatrix<Ta,UnitDiag,RowMajor> tempA = A;
	    TriLDivEq(tempA,B);
	  } else {
	    UpperTriMatrix<Ta,UnitDiag,ColMajor> tempA = A;
	    TriLDivEq(tempA,B);
	  }
	}
      } else {
#ifdef BLAS
	if ( (A.isrm() || A.iscm() ) && (B.isrm() || B.iscm()) ) 
	  BlasTriLDivEq(A,B);
	else 
#endif
	  NonBlasTriLDivEq(A,B);
      }
    }
#ifdef XDEBUG
    Matrix<T> BB = A*B;
    if (Norm(BB-B0) > 0.001*max(RealType(T)(1),Norm(B0))) {
      cerr<<"TriLDivEq: M/Upper\n";
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"Done: B = "<<B<<endl;
      cerr<<"A*B = "<<BB<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta> inline void TriLDivEq(
      const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
#ifdef XDEBUG
    Matrix<T> B0 = B;
#endif

    TMVAssert(A.size() == B.colsize());
    if (B.colsize() > 0 && B.rowsize() > 0) {
      if (B.isconj()) TriLDivEq(A.Conjugate(),B.Conjugate());
      else if (B.rowsize() == 1) TriLDivEq(A,B.col(0));
      else if (B.SameStorageAs(A)) {
	if (A.dt() == NonUnitDiag) {
	  if (A.isrm()) {
	    LowerTriMatrix<Ta,NonUnitDiag,RowMajor> tempA = A;
	    TriLDivEq(tempA,B);
	  } else {
	    LowerTriMatrix<Ta,NonUnitDiag,ColMajor> tempA = A;
	    TriLDivEq(tempA,B);
	  }
	} else {
	  if (A.isrm()) {
	    LowerTriMatrix<Ta,UnitDiag,RowMajor> tempA = A;
	    TriLDivEq(tempA,B);
	  } else {
	    LowerTriMatrix<Ta,UnitDiag,ColMajor> tempA = A;
	    TriLDivEq(tempA,B);
	  }
	}
      } else {
#ifdef BLAS
	if ( (A.isrm() || A.iscm() ) && (B.isrm() || B.iscm()) ) 
	  BlasTriLDivEq(A,B);
	else 
#endif
	  NonBlasTriLDivEq(A,B);
      }
    }
#ifdef XDEBUG
    Matrix<T> BB = A*B;
    if (Norm(BB-B0) > 0.001*max(RealType(T)(1),Norm(B0))) {
      cerr<<"TriLDivEq: M/Lower\n";
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"Done: B = "<<B<<endl;
      cerr<<"A*B = "<<BB<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_TriDiv_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


