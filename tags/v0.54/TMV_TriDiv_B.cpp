
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
  
  //
  // TriLDivEq M
  //

  template <class T, class Ta> inline void RowTriLDivEq(
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
	if (*Aii==Ta(0)) 
	  throw SingularUpperTriMatrix<Ta>(A);
	B.row(i) /= (A.isconj() ? CONJ(*Aii) : *Aii);
      }
    }
  }

  template <class T, class Ta> inline void ColTriLDivEq(
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
	if (*Ajj==Ta(0)) 
	  throw SingularUpperTriMatrix<Ta>(A);
	B.row(j) /= (A.isconj() ? CONJ(*Ajj) : *Ajj);
	B.Rows(0,j) -= A.col(j,0,j) ^ B.row(j);
      }
    } 
  }

  template <class T, class Ta> inline void RowTriLDivEq(
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
	if (*Aii==Ta(0)) 
	  throw SingularLowerTriMatrix<Ta>(A);
	B.row(i) /= (A.isconj() ? CONJ(*Aii) : *Aii);
      }
    }
  }

  template <class T, class Ta> inline void ColTriLDivEq(
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
	if (*Ajj==Ta(0)) 
	  throw SingularLowerTriMatrix<Ta>(A);
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
      return NonBlasTriLDivEq(A11,B1);
    }
  }

#ifdef BLAS
  template <class T, class Ta> inline void BlasTriLDivEq(
      const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
  { NonBlasTriLDivEq(A,B); }
  template <class T, class Ta> inline void BlasTriLDivEq(
      const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
  { NonBlasTriLDivEq(A,B); }
  template <> inline void BlasTriLDivEq(
      const GenUpperTriMatrix<double>& A, const MatrixView<double>& B)
  { 
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.colsize()>0);
    TMVAssert(B.rowsize()>0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(B.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());

    int m=B.iscm()?B.colsize():B.rowsize();
    int n=B.iscm()?B.rowsize():B.colsize();
    double alpha = 1.;
    int lda = A.isrm()?A.stepi():A.stepj();
    int ldb = B.isrm()?B.stepi():B.stepj();

    BLASNAME(dtrsm) (BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
	A.iscm()?BLASCH_UP:BLASCH_LO, 
	A.iscm()==B.iscm()?BLASCH_NT:BLASCH_T,
	A.isunit()?BLASCH_U:BLASCH_NU, 
	BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
	BLASP(B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
  }
  template <> inline void BlasTriLDivEq(
      const GenLowerTriMatrix<double>& A, const MatrixView<double>& B)
  { 
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.colsize()>0);
    TMVAssert(B.rowsize()>0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(B.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());

    int m=B.iscm()?B.colsize():B.rowsize();
    int n=B.iscm()?B.rowsize():B.colsize();
    double alpha = 1.;
    int lda = A.isrm()?A.stepi():A.stepj();
    int ldb = B.isrm()?B.stepi():B.stepj();

    BLASNAME(dtrsm) (BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
	A.iscm()?BLASCH_LO:BLASCH_UP, 
	A.iscm()==B.iscm()?BLASCH_NT:BLASCH_T,
	A.isunit()?BLASCH_U:BLASCH_NU, 
	BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
	BLASP(B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
  }
  template <> inline void BlasTriLDivEq(
      const GenUpperTriMatrix<complex<double> >& A,
      const MatrixView<complex<double> >& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.colsize()>0);
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());

    int m=B.iscm()?B.colsize():B.rowsize();
    int n=B.iscm()?B.rowsize():B.colsize();
    complex<double> alpha = 1.;
    int lda = A.isrm()?A.stepi():A.stepj();
    int ldb = B.isrm()?B.stepi():B.stepj();
    if (A.iscm()==B.iscm() && A.isconj()) {
      B.ConjugateSelf();
      BLASNAME(ztrsm) (BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
	  A.iscm()?BLASCH_UP:BLASCH_LO, BLASCH_NT,
	  A.isunit()?BLASCH_U:BLASCH_NU, 
	  BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
	  BLASP(B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
      B.ConjugateSelf();
    } else {
      BLASNAME(ztrsm) (BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
	  A.iscm()?BLASCH_UP:BLASCH_LO, 
	  A.iscm()==B.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
	  A.isunit()?BLASCH_U:BLASCH_NU, 
	  BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
	  BLASP(B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
    }
  }
  template <> inline void BlasTriLDivEq(
      const GenLowerTriMatrix<complex<double> >& A,
      const MatrixView<complex<double> >& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.colsize()>0);
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());

    int m=B.iscm()?B.colsize():B.rowsize();
    int n=B.iscm()?B.rowsize():B.colsize();
    complex<double> alpha = 1.;
    int lda = A.isrm()?A.stepi():A.stepj();
    int ldb = B.isrm()?B.stepi():B.stepj();
    if (A.iscm()==B.iscm() && A.isconj()) {
      B.ConjugateSelf();
      BLASNAME(ztrsm) (BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
	  A.iscm()?BLASCH_LO:BLASCH_UP, BLASCH_NT,
	  A.isunit()?BLASCH_U:BLASCH_NU, 
	  BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
	  BLASP(B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
      B.ConjugateSelf();
    } else {
      BLASNAME(ztrsm) (BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
	  A.iscm()?BLASCH_LO:BLASCH_UP, 
	  A.iscm()==B.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
	  A.isunit()?BLASCH_U:BLASCH_NU, 
	  BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
	  BLASP(B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
    }
  }
#ifndef NOFLOAT
  template <> inline void BlasTriLDivEq(
      const GenUpperTriMatrix<float>& A, const MatrixView<float>& B)
  { 
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.colsize()>0);
    TMVAssert(B.rowsize()>0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(B.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());

    int m=B.iscm()?B.colsize():B.rowsize();
    int n=B.iscm()?B.rowsize():B.colsize();
    float alpha = 1.;
    int lda = A.isrm()?A.stepi():A.stepj();
    int ldb = B.isrm()?B.stepi():B.stepj();

    BLASNAME(strsm) (BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
	A.iscm()?BLASCH_UP:BLASCH_LO, 
	A.iscm()==B.iscm()?BLASCH_NT:BLASCH_T,
	A.isunit()?BLASCH_U:BLASCH_NU, 
	BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
	BLASP(B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
  }
  template <> inline void BlasTriLDivEq(
      const GenLowerTriMatrix<float>& A, const MatrixView<float>& B)
  { 
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.colsize()>0);
    TMVAssert(B.rowsize()>0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(B.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());

    int m=B.iscm()?B.colsize():B.rowsize();
    int n=B.iscm()?B.rowsize():B.colsize();
    float alpha = 1.;
    int lda = A.isrm()?A.stepi():A.stepj();
    int ldb = B.isrm()?B.stepi():B.stepj();

    BLASNAME(strsm) (BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
	A.iscm()?BLASCH_LO:BLASCH_UP, 
	A.iscm()==B.iscm()?BLASCH_NT:BLASCH_T,
	A.isunit()?BLASCH_U:BLASCH_NU, 
	BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
	BLASP(B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
  }
  template <> inline void BlasTriLDivEq(
      const GenUpperTriMatrix<complex<float> >& A,
      const MatrixView<complex<float> >& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.colsize()>0);
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());

    int m=B.iscm()?B.colsize():B.rowsize();
    int n=B.iscm()?B.rowsize():B.colsize();
    complex<float> alpha = 1.;
    int lda = A.isrm()?A.stepi():A.stepj();
    int ldb = B.isrm()?B.stepi():B.stepj();
    if (A.iscm()==B.iscm() && A.isconj()) {
      B.ConjugateSelf();
      BLASNAME(ctrsm) (BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
	  A.iscm()?BLASCH_UP:BLASCH_LO, BLASCH_NT,
	  A.isunit()?BLASCH_U:BLASCH_NU, 
	  BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
	  BLASP(B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
      B.ConjugateSelf();
    } else {
      BLASNAME(ctrsm) (BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
	  A.iscm()?BLASCH_UP:BLASCH_LO, 
	  A.iscm()==B.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
	  A.isunit()?BLASCH_U:BLASCH_NU, 
	  BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
	  BLASP(B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
    }
  }
  template <> inline void BlasTriLDivEq(
      const GenLowerTriMatrix<complex<float> >& A,
      const MatrixView<complex<float> >& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.colsize()>0);
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());

    int m=B.iscm()?B.colsize():B.rowsize();
    int n=B.iscm()?B.rowsize():B.colsize();
    complex<float> alpha = 1.;
    int lda = A.isrm()?A.stepi():A.stepj();
    int ldb = B.isrm()?B.stepi():B.stepj();
    if (A.iscm()==B.iscm() && A.isconj()) {
      B.ConjugateSelf();
      BLASNAME(ctrsm) (BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
	  A.iscm()?BLASCH_LO:BLASCH_UP, BLASCH_NT,
	  A.isunit()?BLASCH_U:BLASCH_NU, 
	  BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
	  BLASP(B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
      B.ConjugateSelf();
    } else {
      BLASNAME(ctrsm) (BLASCM B.iscm()?BLASCH_L:BLASCH_R, 
	  A.iscm()?BLASCH_LO:BLASCH_UP, 
	  A.iscm()==B.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
	  A.isunit()?BLASCH_U:BLASCH_NU, 
	  BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
	  BLASP(B.ptr()),BLASV(ldb) BLAS1 BLAS1 BLAS1 BLAS1);
    }
  }
#endif // FLOAT
#endif // BLAS

  template <class T, class Ta> void TriLDivEq(
      const GenUpperTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
#ifdef XDEBUG
    Matrix<T> A0 = A;
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
	if (IsComplex(T()) && IsReal(Ta()))
	  BlasTriLDivEq(A,B);
	else if (!((A.isrm()&&A.stepi()>0) || (A.iscm()&&A.stepj()>0) ) ) {
	  if (A.isunit()) {
	    UpperTriMatrix<Ta,UnitDiag,ColMajor> AA(A);
	    TriLDivEq(AA,B);
	  } else {
	    UpperTriMatrix<Ta,NonUnitDiag,ColMajor> AA(A);
	    TriLDivEq(AA,B);
	  }
	} else if (!((B.isrm()&&B.stepi()>0) || (B.iscm()&&B.stepj()>0) ) ) {
	  Matrix<T,ColMajor> BB(B);
	  BlasTriLDivEq(A,BB.View());
	  B = BB;
	} else 
	  BlasTriLDivEq(A,B);
#else
	NonBlasTriLDivEq(A,B);
#endif
      }
    }
#ifdef XDEBUG
    Matrix<T> BB = A0*B;
    if (Norm(BB-B0) > 0.001*Norm(A0)*Norm(B0)) {
      cerr<<"TriLDivEq: M/Upper\n";
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"Done: B = "<<B<<endl;
      cerr<<"A*B = "<<BB<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta> void TriLDivEq(
      const GenLowerTriMatrix<Ta>& A, const MatrixView<T>& B)
  {
#ifdef XDEBUG
    Matrix<T> A0 = A;
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
	if (IsComplex(T()) && IsReal(Ta()))
	  BlasTriLDivEq(A,B);
	else if (!((A.isrm()&&A.stepi()>0) || (A.iscm()&&A.stepj()>0) ) ) {
	  if (A.isunit()) {
	    LowerTriMatrix<Ta,UnitDiag,ColMajor> AA(A);
	    TriLDivEq(AA,B);
	  } else {
	    LowerTriMatrix<Ta,NonUnitDiag,ColMajor> AA(A);
	    TriLDivEq(AA,B);
	  }
	} else if (!((B.isrm()&&B.stepi()>0) || (B.iscm()&&B.stepj()>0) ) ) {
	  Matrix<T,ColMajor> BB(B);
	  BlasTriLDivEq(A,BB.View());
	  B = BB;
	} else 
	  BlasTriLDivEq(A,B);
#else
	NonBlasTriLDivEq(A,B);
#endif
      }
    }
#ifdef XDEBUG
    Matrix<T> BB = A0*B;
    if (Norm(BB-B0) > 0.001*Norm(A0)*Norm(B0)) {
      cerr<<"TriLDivEq: M/Lower\n";
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
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


