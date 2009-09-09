
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
  // TriLDivEq V
  //

  template <bool rm, bool ca, bool ua, class T, class Ta> void DoRowTriLDivEq(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& b)
  {
    // Solve A x = y  where A is an upper triangle matrix
    //cerr<<"Row Upper\n";
    TMVAssert(b.step()==1);
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size() > 0);
    TMVAssert(b.ct() == NonConj);
    TMVAssert(rm == A.isrm());
    TMVAssert(ca == A.isconj());
    TMVAssert(ua == A.isunit());

    const size_t N = A.size();

    const int sj = (rm?1:A.stepj());
    const int ds = A.stepi()+sj;
    const Ta* Aii = A.cptr() + (ua ? N-2 : N-1)*ds;
    T* bi = b.ptr() + (ua ? N-2 : N-1);

    if (!ua) {
      if (*Aii==Ta(0)) tmv_error("Singular Matrix found in UpperTriLDivEq");
      *bi /= (ca ? CONJ(*Aii) : *Aii);
      Aii -= ds;
      --bi;
    }
    if (N==1) return;

    for(size_t i=N-1,len=1; i>0; --i,++len,Aii-=ds,--bi) {
      // Actual row being done is i-1, not i

      // *bi -= A.row(i,i+1,N) * b.SubVector(i+1,N);
      const T* bj = bi+1;
      const Ta* Aij = Aii + sj;
      for(size_t j=len;j>0;--j,++bj,(rm?++Aij:Aij+=sj))
	*bi -= (*bj) * (ca ? CONJ(*Aij) : *Aij);

      if (!ua) {
	if (*Aii==Ta(0)) 
	  tmv_error("Singular Matrix found in UpperTriLDivEq");
	*bi /= (ca ? CONJ(*Aii) : *Aii);
      }
    }
  }

  template <bool rm, class T, class Ta> void RowTriLDivEq(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& b)
  {
    if (A.isconj())
      if (A.isunit())
	DoRowTriLDivEq<rm,true,true>(A,b);
      else
	DoRowTriLDivEq<rm,true,false>(A,b);
    else
      if (A.isunit())
	DoRowTriLDivEq<rm,false,true>(A,b);
      else
	DoRowTriLDivEq<rm,false,false>(A,b);
  }

  template <bool cm, bool ca, bool ua, class T, class Ta> void DoColTriLDivEq(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& b)
  {
    //cerr<<"colmajor upper\n";
    // Solve A x = y  where A is an upper triangle matrix
    TMVAssert(b.step()==1);
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size() > 0);
    TMVAssert(b.ct() == NonConj);
    TMVAssert(cm == A.iscm());
    TMVAssert(ca == A.isconj());
    TMVAssert(ua == A.isunit());

    const size_t N = A.size();

    const int si = (cm ? 1 : A.stepi());
    const int sj = A.stepj();
    const int ds = si+sj;
    const Ta* A0j = A.cptr()+(N-1)*sj;
    const Ta* Ajj = (ua ? 0 : A0j+(N-1)*si); // if unit, this isn't used.
    T*const b0 = b.ptr();
    T* bj = b0 + N-1;

    for(size_t j=N-1; j>0; --j,--bj,A0j-=sj) {
      if (*bj != T(0)) {
	if (!ua) {
	  if (*Ajj==Ta(0)) tmv_error("Singular Matrix found in UpperTriLDivEq");
	  *bj /= (ca ? CONJ(*Ajj) : *Ajj);
	  Ajj-=ds;
	}

	// b.SubVector(0,j) -= *bj * A.col(j,0,j);
	T* bi = b0;
	const Ta* Aij = A0j;
	for(size_t i=j;i>0;--i,++bi,(cm?++Aij:Aij+=si))
	  *bi -= *bj * (ca ? CONJ(*Aij) : *Aij);
      }
      else if (!ua) Ajj -= ds;
    }
    if (!ua && *bj != T(0)) {
      if (*Ajj==Ta(0)) tmv_error("Singular Matrix found in UpperTriLDivEq");
      *bj /= (ca?CONJ(*Ajj):*Ajj);
    } 
  }

  template <bool cm, class T, class Ta> void ColTriLDivEq(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& b)
  {
    if (A.isconj())
      if (A.isunit())
	DoColTriLDivEq<cm,true,true>(A,b);
      else
	DoColTriLDivEq<cm,true,false>(A,b);
    else
      if (A.isunit())
	DoColTriLDivEq<cm,false,true>(A,b);
      else
	DoColTriLDivEq<cm,false,false>(A,b);
  }

  template <bool rm, bool ca, bool ua, class T, class Ta> void DoRowTriLDivEq(
	const GenLowerTriMatrix<Ta>& A, const VectorView<T>& b)
  {
    // Solve A x = y  where A is a lower triangle matrix
    TMVAssert(b.step()==1);
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size() > 0);
    TMVAssert(b.ct() == NonConj);
    TMVAssert(rm == A.isrm());
    TMVAssert(ca == A.isconj());
    TMVAssert(ua == A.isunit());

    const size_t N = A.size();

    const int sj = (rm ? 1 : A.stepj());
    const int si = A.stepi();

    const Ta* Ai0 = A.cptr();
    T* b0 = b.ptr();

    if (!ua) {
      if (*Ai0==Ta(0)) tmv_error("Singular Matrix found in LowerTriLDivEq");
      *b0 /= (ca ? CONJ(*Ai0) : *Ai0);
    }

    T* bi = b0+1;
    Ai0 += si;
    for(size_t i=1,len=1;i<N;++i,++len,++bi,Ai0+=si) {
      // *bi -= A.row(i,0,i) * b.SubVector(0,i);
      const Ta* Aij = Ai0;
      const T* bj = b0;
      for(size_t j=len;j>0;--j,++bj,(rm?++Aij:Aij+=sj))
	*bi -= (*bj) * (ca ? CONJ(*Aij) : *Aij);
      if (!ua) {
	// Aij is Aii after the above for loop
	if (*Aij==Ta(0)) 
	  tmv_error("Singular Matrix found in LowerTriLDivEq");
	*bi /= (ca ? CONJ(*Aij) : *Aij);
      }
    }
  }

  template <bool rm, class T, class Ta> void RowTriLDivEq(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& b)
  {
    if (A.isconj())
      if (A.isunit())
	DoRowTriLDivEq<rm,true,true>(A,b);
      else
	DoRowTriLDivEq<rm,true,false>(A,b);
    else
      if (A.isunit())
	DoRowTriLDivEq<rm,false,true>(A,b);
      else
	DoRowTriLDivEq<rm,false,false>(A,b);
  }

  template <bool cm, bool ca, bool ua, class T, class Ta> void DoColTriLDivEq(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& b)
  {
    // Solve A x = y  where A is a lower triangle matrix
    TMVAssert(b.step()==1);
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size() > 0);
    TMVAssert(b.ct() == NonConj);
    TMVAssert(cm == A.iscm());
    TMVAssert(ca == A.isconj());
    TMVAssert(ua == A.isunit());

    const size_t N = A.size();

    const int si = (cm ? 1 : A.stepi());
    const int ds = A.stepj()+si;
    const Ta* Ajj = A.cptr();
    T* bj = b.ptr();

    for(size_t j=0,len=N-1;len>0;++j,--len,++bj,Ajj+=ds) if (*bj != T(0)) {
      if (!ua) {
	if (*Ajj==Ta(0)) tmv_error("Singular Matrix found in LowerTriLDivEq");
	*bj /= (ca ? CONJ(*Ajj) : *Ajj);
      }
      // b.SubVecotr(j+1,N) -= *bj * A.col(j,j+1,N);
      T* bi = bj+1;
      const Ta* Aij = Ajj+si;
      for(size_t i=len;i>0;--i,++bi,(cm?++Aij:Aij+=si))
	*bi -= *bj * (ca ? CONJ(*Aij) : *Aij);
    }
    if (!ua && *bj != T(0)) {
      if (*Ajj==Ta(0)) tmv_error("Singular Matrix found in LowerTriLDivEq");
      *bj /= (ca ? CONJ(*Ajj) : *Ajj);
    } 
  }

  template <bool cm, class T, class Ta> void ColTriLDivEq(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& b)
  {
    if (A.isconj())
      if (A.isunit())
	DoColTriLDivEq<cm,true,true>(A,b);
      else
	DoColTriLDivEq<cm,true,false>(A,b);
    else
      if (A.isunit())
	DoColTriLDivEq<cm,false,true>(A,b);
      else
	DoColTriLDivEq<cm,false,false>(A,b);
  }

  template <class T, class Ta> void DoTriLDivEq(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& b)
  {
    if (A.isrm()) RowTriLDivEq<true>(A,b);
    else if (A.iscm()) ColTriLDivEq<true>(A,b);
    else RowTriLDivEq<false>(A,b); 
  }

  template <class T, class Ta> void DoTriLDivEq(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& b)
  {
    if (A.isrm()) RowTriLDivEq<true>(A,b);
    else if (A.iscm()) ColTriLDivEq<true>(A,b);
    else RowTriLDivEq<false>(A,b); 
  }

  template <class T, class Ta> void NonBlasTriLDivEq(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& b)
  {
    //cerr<<"Upper LDivEq vect\n";
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(b.ct() == NonConj);

    if (b.step() == 1) {
      size_t i2 = b.size();
      for(const T* b2 = b.cptr()+i2-1; i2>0 && *b2==T(0); --i2,--b2);
      if (i2==0) return;
      else if (i2 == b.size())
	DoTriLDivEq(A,b);
      else
	DoTriLDivEq(A.SubTriMatrix(0,i2),b.SubVector(0,i2));
    } else {
      Vector<T> bb = b;
      NonBlasTriLDivEq(A,bb.View());
      b = bb;
    }
  }

  template <class T, class Ta> void NonBlasTriLDivEq(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& b)
  {
    //cerr<<"Lower LDivEq vect\n";
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(b.ct() == NonConj);

    if (b.step() == 1) {
      const size_t N = b.size();
      size_t i1 = 0;
      for(const T* b1 = b.cptr(); i1<N && *b1==T(0); ++i1,++b1);
      if (i1==N) return;
      else if (i1 == 0)
	DoTriLDivEq(A,b);
      else
	DoTriLDivEq(A.SubTriMatrix(i1,N),b.SubVector(i1,N));
    } else {
      Vector<T> bb = b;
      NonBlasTriLDivEq(A,bb.View());
      b = bb;
    }
  }

#ifdef BLAS
  template <class T, class Ta> inline void BlasTriLDivEq(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& b)
  { NonBlasTriLDivEq(A,b); }
  template <class T, class Ta> inline void BlasTriLDivEq(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& b)
  { NonBlasTriLDivEq(A,b); }
  inline void BlasTriLDivEq(const GenUpperTriMatrix<double>& A,
      const VectorView<double>& b)
  {
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(b.ct() == NonConj);
    cblas_dtrsv(A.isrm() ? CblasRowMajor : CblasColMajor,
	CblasUpper, CblasNoTrans,
	A.isunit() ? CblasUnit : CblasNonUnit,
	A.size(),A.cptr(),A.isrm() ? A.stepi() : A.stepj(),
	b.ptr(),b.step());
  }
  inline void BlasTriLDivEq(const GenLowerTriMatrix<double>& A,
      const VectorView<double>& b)
  {
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(b.ct() == NonConj);
    cblas_dtrsv(A.isrm() ? CblasRowMajor : CblasColMajor,
	CblasLower, CblasNoTrans,
	A.isunit() ? CblasUnit : CblasNonUnit,
	A.size(),A.cptr(),A.isrm() ? A.stepi() : A.stepj(),
	b.ptr(),b.step());
  }
  inline void BlasTriLDivEq(const GenUpperTriMatrix<complex<double> >& A,
      const VectorView<complex<double> >& b)
  {
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(b.ct() == NonConj);
    CBLAS_ORDER order = (A.iscm() == A.isconj()) ?
      CblasRowMajor : CblasColMajor;
    CBLAS_UPLO uplo = A.isconj() ?  CblasLower : CblasUpper;
    CBLAS_TRANSPOSE tran = A.isconj() ?  CblasConjTrans : CblasNoTrans;
    cblas_ztrsv(order,uplo,tran, A.isunit() ? CblasUnit : CblasNonUnit,
	A.size(),A.cptr(),A.isrm() ? A.stepi() : A.stepj(),
	b.ptr(),b.step());
  }
  inline void BlasTriLDivEq(const GenLowerTriMatrix<complex<double> >& A,
      const VectorView<complex<double> >& b)
  {
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(b.ct() == NonConj);
    CBLAS_ORDER order = (A.iscm() == A.isconj()) ?
      CblasRowMajor : CblasColMajor;
    CBLAS_UPLO uplo = A.isconj() ?  CblasUpper : CblasLower;
    CBLAS_TRANSPOSE tran = A.isconj() ?  CblasConjTrans : CblasNoTrans;
    cblas_ztrsv(order,uplo,tran, A.isunit() ? CblasUnit : CblasNonUnit,
	A.size(),A.cptr(),A.isrm() ? A.stepi() : A.stepj(),
	b.ptr(),b.step());
  }
#ifndef NOFLOAT
  inline void BlasTriLDivEq(const GenUpperTriMatrix<float>& A,
      const VectorView<float>& b)
  {
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(b.ct() == NonConj);
    cblas_strsv(A.isrm() ? CblasRowMajor : CblasColMajor,
	CblasUpper, CblasNoTrans,
	A.isunit() ? CblasUnit : CblasNonUnit,
	A.size(),A.cptr(),A.isrm() ? A.stepi() : A.stepj(),
	b.ptr(),b.step());
  }
  inline void BlasTriLDivEq(const GenLowerTriMatrix<float>& A,
      const VectorView<float>& b)
  {
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(b.ct() == NonConj);
    cblas_strsv(A.isrm() ? CblasRowMajor : CblasColMajor,
	CblasLower, CblasNoTrans,
	A.isunit() ? CblasUnit : CblasNonUnit,
	A.size(),A.cptr(),A.isrm() ? A.stepi() : A.stepj(),
	b.ptr(),b.step());
  }
  inline void BlasTriLDivEq(const GenUpperTriMatrix<complex<float> >& A,
      const VectorView<complex<float> >& b)
  {
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(b.ct() == NonConj);
    CBLAS_ORDER order = (A.iscm() == A.isconj()) ?
      CblasRowMajor : CblasColMajor;
    CBLAS_UPLO uplo = A.isconj() ?  CblasLower : CblasUpper;
    CBLAS_TRANSPOSE tran = A.isconj() ?  CblasConjTrans : CblasNoTrans;
    cblas_ctrsv(order,uplo,tran, A.isunit() ? CblasUnit : CblasNonUnit,
	A.size(),A.cptr(),A.isrm() ? A.stepi() : A.stepj(),
	b.ptr(),b.step());
  }
  inline void BlasTriLDivEq(const GenLowerTriMatrix<complex<float> >& A,
      const VectorView<complex<float> >& b)
  {
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(b.ct() == NonConj);
    CBLAS_ORDER order = (A.iscm() == A.isconj()) ?
      CblasRowMajor : CblasColMajor;
    CBLAS_UPLO uplo = A.isconj() ?  CblasUpper : CblasLower;
    CBLAS_TRANSPOSE tran = A.isconj() ?  CblasConjTrans : CblasNoTrans;
    cblas_ctrsv(order,uplo,tran, A.isunit() ? CblasUnit : CblasNonUnit,
	A.size(),A.cptr(),A.isrm() ? A.stepi() : A.stepj(),
	b.ptr(),b.step());
  }
#endif
#endif

  template <class T, class Ta> inline void TriLDivEq(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& b)
  {
    TMVAssert(b.size() == A.size());
#ifdef XDEBUG
    Vector<T> b0 = b;
#endif
    if (b.size() > 0) {
      if (b.isconj()) TriLDivEq(A.Conjugate(),b.Conjugate());
      else 
#ifdef BLAS
	if ( A.isrm() || A.iscm() )
	  BlasTriLDivEq(A,b);
	else 
#endif
	  NonBlasTriLDivEq(A,b);
    }
#ifdef XDEBUG
    Vector<T> b2 = A*b;
    if (Norm(b2-b0) > 0.001*max(RealType(T)(1),Norm(b0))) {
      cerr<<"TriLDivEq: v/Upper\n";
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"b = "<<Type(b)<<"  "<<b0<<endl;
      cerr<<"Done: b = "<<b<<endl;
      cerr<<"A*b = "<<b2<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta> inline void TriLDivEq(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& b)
  {
    TMVAssert(b.size() == A.size());
#ifdef XDEBUG
    Vector<T> b0 = b;
#endif
    if (b.size() > 0) {
      if (b.isconj()) TriLDivEq(A.Conjugate(),b.Conjugate());
      else 
#ifdef BLAS
	if ( A.isrm() || A.iscm() )
	  BlasTriLDivEq(A,b);
	else 
#endif
	  NonBlasTriLDivEq(A,b);
    }
#ifdef XDEBUG
    Vector<T> b2 = A*b;
    if (Norm(b2-b0) > 0.001*max(RealType(T)(1),Norm(b0))) {
      cerr<<"TriLDivEq: v/Lower\n";
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"b = "<<Type(b)<<"  "<<b0<<endl;
      cerr<<"Done: b = "<<b<<endl;
      cerr<<"A*b = "<<b2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_TriDiv_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


