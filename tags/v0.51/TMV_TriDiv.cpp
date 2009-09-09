
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
  
  template <class T> bool UpperTriDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenUpperTriMatrix<T>* tm = dynamic_cast<const GenUpperTriMatrix<T>*>(&m);
    TMVAssert(tm);
    if (fout) {
      *fout << "M = "<<tmv::Type(*tm)<<"  "<<*tm<<endl;
      *fout << "T = "<<*itsm<<endl;
    }
    RealType(T) nm = Norm(*itsm-*tm);
    nm /= Norm(*itsm);
    if (fout) {
      *fout << "Norm(M-T)/Norm(T) = "<<nm<<endl;
    }
    return nm < Epsilon<T>();
  }

  template <class T> bool LowerTriDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenLowerTriMatrix<T>* tm = dynamic_cast<const GenLowerTriMatrix<T>*>(&m);
    TMVAssert(tm);
    if (fout) {
      *fout << "M = "<<tmv::Type(*tm)<<"  "<<*tm<<endl;
      *fout << "T = "<<*itsm<<endl;
    }
    RealType(T) nm = Norm(*itsm-*tm);
    nm /= Norm(*itsm);
    if (fout) {
      *fout << "Norm(M-T)/Norm(T) = "<<nm<<endl;
    }
    return nm < Epsilon<T>();
  }

  template <class T> T UpperTriDiv<T>::Det() const
  {
    if (!donedet) {
      if (!itsm->isunit()) det *= DiagMatrixViewOf(itsm->diag()).Det();
      donedet = true;
    }
    return det;  
  }                  

  template <class T> T LowerTriDiv<T>::Det() const
  {
    if (!donedet) {
      if (!itsm->isunit()) det *= DiagMatrixViewOf(itsm->diag()).Det();
      donedet = true;
    }
    return det;  
  }                  

  template <class T> void UpperTriDiv<T>::TInverse(
      const UpperTriMatrixView<T>& minv) const
  {
    TMVAssert(minv.size() == itsm->size());
    TMVAssert(!minv.isunit() || itsm->isunit());
    minv.SetToIdentity();
    LDivEq(minv);
  }

  template <class T> void LowerTriDiv<T>::TInverse(
      const LowerTriMatrixView<T>& minv) const
  {
    TMVAssert(minv.size() == itsm->size());
    TMVAssert(!minv.isunit() || itsm->isunit());
    minv.SetToIdentity();
    LDivEq(minv);
  }

  template <class T> void UpperTriDiv<T>::InverseATA(
      const MatrixView<T>& minv) const
  {
    TMVAssert(minv.colsize() == itsm->size());
    TMVAssert(minv.rowsize() == itsm->size());
    UpperTriMatrix<T,NonUnitDiag,ColMajor> temp(itsm->size());
    TInverse(temp.View());
    minv = temp*Adjoint(temp);
  }

  template <class T> void LowerTriDiv<T>::InverseATA(
      const MatrixView<T>& minv) const
  {
    TMVAssert(minv.colsize() == itsm->size());
    TMVAssert(minv.rowsize() == itsm->size());
    LowerTriMatrix<T,NonUnitDiag,ColMajor> temp(itsm->size());
    TInverse(temp.View());
    minv = temp*Adjoint(temp);
  }

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
#endif
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
	if (*Aii==Ta(0)) tmv_error("Singular Matrix found in UpperTriLDivEq");
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
	  tmv_error("Singular Matrix found in UpperTriLDivEq");
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
	if (*Aii==Ta(0)) tmv_error("Singular Matrix found in LowerTriLDivEq");
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
	  tmv_error("Singular Matrix found in LowerTriLDivEq");
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

  template <class T, class Ta> inline void TriLDivEq(
      const GenUpperTriMatrix<Ta>& A, const UpperTriMatrixView<T>& B)
    // B = A^-1 * B
    // where A is a triangle matrix
  {
#ifdef XDEBUG
    UpperTriMatrix<T> B0 = B;
#endif
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || A.isunit());

    if (B.isconj()) TriLDivEq(A.Conjugate(),B.Conjugate());
    else DoTriLDivEq(A,B);

#ifdef XDEBUG
    UpperTriMatrix<T> BB = A*B;
    if (Norm(BB-B0) > 0.001*max(RealType(T)(1),Norm(B0))) {
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

  template <class T, class Ta> inline void TriLDivEq(
      const GenLowerTriMatrix<Ta>& A, const LowerTriMatrixView<T>& B)
    // B = A^-1 * B
    // where A is a triangle matrix
  {
#ifdef XDEBUG
    LowerTriMatrix<T> B0 = B;
#endif

    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || A.isunit());

    if (B.isconj()) TriLDivEq(A.Conjugate(),B.Conjugate());
    else if (B.size() > 0) DoTriLDivEq(A,B);

#ifdef XDEBUG
    LowerTriMatrix<T> BB = A*B;
    if (Norm(BB-B0) > 0.001*max(RealType(T)(1),Norm(B0))) {
      cerr<<"TriLDivEq: Lower/Lower\n";
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"Done: B = "<<B<<endl;
      cerr<<"A*B = "<<BB<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_TriDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


