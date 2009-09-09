
#include "TMV_VectorArith_Inline.h"
#include "TMV_Tri.h"
#include "TMV_Diag.h"

//#define XDEBUG

namespace tmv {

  // MJ: Compare with ATLAS - try to speed up
  
  template <class T> bool UpperTriDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenUpperTriMatrix<T>* tm = dynamic_cast<const GenUpperTriMatrix<T>*>(&m);
    TMVAssert(tm);
    if (fout) {
      *fout << "M = "<<*tm<<endl;
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
      *fout << "M = "<<*tm<<endl;
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

  template <class T1, class T> void RowMajorTriLDivEq(
	const GenUpperTriMatrix<T1>& A, const VectorView<T>& b)
  {
    //cerr<<"rowmajor upper\n";
    // Solve A x = y  where A is an upper triangle matrix
    TMVAssert(A.isrm());
    TMVAssert(b.step()==1);
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(b.ct() == NonConj);

    size_t N = A.size();
    VIt<T,Unit,NonConj> bi = b.begin()+N-1;
    for(;N>0 && *bi==T(0);--N,--bi);
    if (N==0) return;

    const int ds = A.stepi()+1;
    const T1* Aptr = A.cptr()+(N-2)*ds+1;

    if (A.isunit()) {
      if (N==1) return;
      for(int i=N-2,len=1; i>=0; --i,++len,Aptr-=ds) {
	// *bi -= A.row(i,i+1,N) * b.SubVector(i+1,N);
	T temp;
	if (A.isconj())
	  DoMultVV(temp,CVIt<T1,Unit,Conj>(Aptr,1),
	      CVIt<T,Unit,NonConj>(bi),len); // here, bi = b(i+1)
	else
	  DoMultVV(temp,CVIt<T1,Unit,NonConj>(Aptr,1),
	      CVIt<T,Unit,NonConj>(bi),len); // here, bi = b(i+1)
	*(--bi) -= temp;
      }
    } else {
      CVIter<T1> Aii = A.diag().begin()+N-1;
      if (*Aii==T1(0)) tmv_error("Singular Matrix found in UpperTriLDivEq");
      *bi /= *Aii;
      for(int i=N-2,len=1; i>=0; --i,++len,Aptr-=ds) {
	// *bi -= A.row(i,i+1,N) * b.SubVector(i+1,N);
	T temp;
	if (A.isconj())
	  DoMultVV(temp,CVIt<T1,Unit,Conj>(Aptr,1),
	      CVIt<T,Unit,NonConj>(bi),len); // here, bi = b(i+1)
	else
	  DoMultVV(temp,CVIt<T1,Unit,NonConj>(Aptr,1),
	      CVIt<T,Unit,NonConj>(bi),len); // here, bi = b(i+1)
	*(--bi) -= temp;
	if (*(--Aii)==T1(0)) 
	  tmv_error("Singular Matrix found in UpperTriLDivEq");
	*bi /= *Aii;
      }
    }
  }

  template <class T1, class T> void RowTriLDivEq(
      const GenUpperTriMatrix<T1>& A, const VectorView<T>& b)
  {
    //cerr<<"row upper\n";
    // Solve A x = y  where A is an upper triangle matrix
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(b.ct() == NonConj);

    size_t N = A.size();
    VIt<T,Step,NonConj> bi = b.begin()+N-1;
    for(;N>0 && *bi==T(0);--N,--bi);
    if (N==0) return;

    if (A.isunit()) {
      for(int i=N-1; i>=0; --i,--bi) 
	*bi -= A.row(i,i+1,N) * b.SubVector(i+1,N);
    } else {
      CVIter<T1> Aii = A.diag().begin()+N-1;
      for(int i=N-1; i>=0; --i,--bi,--Aii) {
	*bi -= A.row(i,i+1,N) * b.SubVector(i+1,N);
	if (*Aii==T1(0)) tmv_error("Singular Matrix found in UpperTriLDivEq");
	*bi /= *Aii;
      }
    }
  }

  template <class T1, class T> void ColMajorTriLDivEq(
      const GenUpperTriMatrix<T1>& A, const VectorView<T>& b)
  {
    //cerr<<"colmajor upper\n";
    // Solve A x = y  where A is an upper triangle matrix
    TMVAssert(A.iscm());
    TMVAssert(b.step()==1);
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(b.ct() == NonConj);

    size_t N = A.size();
    VIt<T,Unit,NonConj> bj = b.begin()+N-1;
    for(;N>0 && *bj==T(0);--N,--bj);
    if (N==0) return;
    const T1* Aptr = A.cptr()+(N-1)*A.stepj();

    if (A.isunit()) {
      for(int j=N-1; j>0; --j,--bj,Aptr-=A.stepj()) if (*bj != T(0))
	// b.SubVector(0,j) -= *bj * A.col(j,0,j);
	if (A.isconj())
	  DoAddVV(-*bj,CVIt<T1,Unit,Conj>(Aptr,1),
	      VIt<T,Unit,NonConj>(b.begin()),j);
	else
	  DoAddVV(-*bj,CVIt<T1,Unit,NonConj>(Aptr,1),
	      VIt<T,Unit,NonConj>(b.begin()),j);
    } else {
      CVIter<T1> Ajj = A.diag().begin()+N-1;
      for(int j=N-1; j>0; --j,--bj,--Ajj,Aptr-=A.stepj()) if (*bj != T(0)) {
	if (*Ajj==T1(0)) tmv_error("Singular Matrix found in UpperTriLDivEq");
	*bj /= *Ajj;
	// b.SubVector(0,j) -= *bj * A.col(j,0,j);
	if (A.isconj())
	  DoAddVV(-*bj,CVIt<T1,Unit,Conj>(Aptr,1),
	      VIt<T,Unit,NonConj>(b.begin()),j);
	else
	  DoAddVV(-*bj,CVIt<T1,Unit,NonConj>(Aptr,1),
	      VIt<T,Unit,NonConj>(b.begin()),j);
      }
      if (*bj != T(0)) {
	if (*Ajj==T1(0)) tmv_error("Singular Matrix found in UpperTriLDivEq");
	*bj /= *Ajj;
      }
    } 
  }

  template <class T1, class T> void RowMajorTriLDivEq(
	const GenLowerTriMatrix<T1>& A, const VectorView<T>& b)
  {
    // Solve A x = y  where A is a lower triangle matrix
    TMVAssert(A.isrm());
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(b.ct() == NonConj);

    const size_t N = A.size();
    size_t i1=0; // all b(i<i1) == 0
    VIt<T,Unit,NonConj> bi1 = b.begin();
    for(;i1 < N && *bi1==T(0);++i1,++bi1);
    if (i1==N) return;

    const T1* Aptr = A.cptr()+i1*(A.stepi()+1)+A.stepi();

    VIt<T,Unit,NonConj> bi = bi1+1;
    if (A.isunit()) {
      for(size_t i=i1+1,len=1;i<N;++i,++len,++bi,Aptr+=A.stepi()) {
	// *bi -= A.row(i,i1,i) * b.SubVector(i1,i);
	T temp;
	if (A.isconj()) 
	  DoMultVV(temp,CVIt<T1,Unit,Conj>(Aptr,1),
	      CVIt<T,Unit,NonConj>(bi1),len);
	else
	  DoMultVV(temp,CVIt<T1,Unit,NonConj>(Aptr,1),
	      CVIt<T,Unit,NonConj>(bi1),len);
	*bi -= temp;
      }
    } else {
      CVIter<T1> Aii = A.diag().begin()+i1;
      if (*Aii==T1(0)) tmv_error("Singular Matrix found in LowerTriLDivEq");
      *bi1 /= *Aii;
      for(size_t i=i1+1,len=1;i<N;++i,++len,++bi,Aptr+=A.stepi()) {
	// *bi -= A.row(i,i1,i) * b.SubVector(i1,i);
	T temp;
	if (A.isconj())
	  DoMultVV(temp,CVIt<T1,Unit,Conj>(Aptr,1),
	      CVIt<T,Unit,NonConj>(bi1),len);
	else
	  DoMultVV(temp,CVIt<T1,Unit,NonConj>(Aptr,1),
	      CVIt<T,Unit,NonConj>(bi1),len);
	*bi -= temp;
	if (*(++Aii)==T1(0)) 
	  tmv_error("Singular Matrix found in LowerTriLDivEq");
	*bi /= *Aii;
      }
    }
  }

  template <class T1, class T> void RowTriLDivEq(
      const GenLowerTriMatrix<T1>& A, const VectorView<T>& b)
  {
    // Solve A x = y  where A is a lower triangle matrix
    const size_t N = A.size();
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(b.ct() == NonConj);

    size_t i1=0; // all b(i<i1) == 0
    VIt<T,Step,NonConj> bi = b.begin();
    for(;i1<N && *bi==T(0);++i1,++bi);
    if (i1==N) return;

    if (A.isunit()) {
      for(size_t i=i1;i<N;++i,++bi) 
	*bi -= A.row(i,i1,i) * b.SubVector(i1,i);
    } else {
      CVIter<T1> Aii = A.diag().begin()+i1;
      for(size_t i=i1;i<N;++i,++Aii,++bi) {
	*bi -= A.row(i,i1,i) * b.SubVector(i1,i);
	if (*Aii==T1(0)) tmv_error("Singular Matrix found in LowerTriLDivEq");
	*bi /= *Aii;
      }
    }
  }

  template <class T1, class T> void ColMajorTriLDivEq(
	const GenLowerTriMatrix<T1>& A, const VectorView<T>& b)
  {
    // Solve A x = y  where A is a lower triangle matrix
    TMVAssert(A.iscm());
    TMVAssert(b.step()==1);
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(b.ct() == NonConj);

    const size_t N = A.size();
    size_t i1=0; // all b(i<i1) == 0
    VIt<T,Unit,NonConj> bj = b.begin();
    for(;i1<N && *bj==T(0);++i1,++bj);
    if (i1==N) return;

    const int ds = A.stepj()+1;
    const T1* Aptr = A.cptr()+i1*ds+1;

    if (A.isunit()) {
      for(size_t j=i1,len=N-j-1;len>0;++j,--len,++bj,Aptr+=ds) {
	// b.SubVecotr(j+1,N) -= *bj * A.col(j,j+1,N);
	if (*bj != T(0)) {
	  if (A.isconj())
	    DoAddVV(-*bj,CVIt<T1,Unit,Conj>(Aptr,1),bj+1,len);
	  else
	    DoAddVV(-*bj,CVIt<T1,Unit,NonConj>(Aptr,1),bj+1,len);
	}
      }
    } else {
      CVIter<T1> Ajj = A.diag().begin()+i1;
      for(size_t j=i1,len=N-j-1;len>0;++j,--len,++Ajj,++bj,Aptr+=ds) {
	if (*bj != T(0)) {
	  if (*Ajj==T1(0)) tmv_error("Singular Matrix found in LowerTriLDivEq");
	  *bj /= *Ajj;
	  // b.SubVecotr(j+1,N) -= *bj * A.col(j,j+1,N);
	  if (A.isconj())
	    DoAddVV(-*bj,CVIt<T1,Unit,Conj>(Aptr,1),bj+1,len);
	  else
	    DoAddVV(-*bj,CVIt<T1,Unit,NonConj>(Aptr,1),bj+1,len);
	}
      }
      if (*bj != T(0)) {
	if (*Ajj==T1(0)) tmv_error("Singular Matrix found in LowerTriLDivEq");
	*bj /= *Ajj;
      } 
    }
  }

  template <class T1, class T> void NonBlasTriLDivEq(
      const GenUpperTriMatrix<T1>& A, const VectorView<T>& b)
  {
    //cerr<<"Upper LDivEq vect\n";
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(b.ct() == NonConj);

    if (b.step() == 1) {
      if (A.iscm()) ColMajorTriLDivEq(A,b);
      else if (A.isrm()) RowMajorTriLDivEq(A,b);
      else RowTriLDivEq(A,b); 
    }
    else RowTriLDivEq(A,b); 
  }

  template <class T1, class T> void NonBlasTriLDivEq(
      const GenLowerTriMatrix<T1>& A, const VectorView<T>& b)
  {
    //cerr<<"Lower LDivEq vect\n";
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(b.ct() == NonConj);

    if (b.step() == 1) {
      if (A.iscm()) ColMajorTriLDivEq(A,b);
      else if (A.isrm()) RowMajorTriLDivEq(A,b);
      else RowTriLDivEq(A,b); 
    }
    else RowTriLDivEq(A,b);
  }

#ifdef BLAS
  template <class T1, class T> inline void BlasTriLDivEq(
      const GenUpperTriMatrix<T1>& A, const VectorView<T>& b)
  { NonBlasTriLDivEq(A,b); }
  template <class T1, class T> inline void BlasTriLDivEq(
      const GenLowerTriMatrix<T1>& A, const VectorView<T>& b)
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
  inline void BlasTriLDivEq(const GenUpperTriMatrix<double>& A,
      const VectorView<complex<double> >& b)
  {
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(b.ct() == NonConj);
    BlasTriLDivEq(A,b.Real());
    BlasTriLDivEq(A,b.Imag());
  }
  inline void BlasTriLDivEq(const GenLowerTriMatrix<double>& A,
      const VectorView<complex<double> >& b)
  {
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(b.ct() == NonConj);
    BlasTriLDivEq(A,b.Real());
    BlasTriLDivEq(A,b.Imag());
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
  inline void BlasTriLDivEq(const GenUpperTriMatrix<float>& A,
      const VectorView<complex<float> >& b)
  {
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(b.ct() == NonConj);
    BlasTriLDivEq(A,b.Real());
    BlasTriLDivEq(A,b.Imag());
  }
  inline void BlasTriLDivEq(const GenLowerTriMatrix<float>& A,
      const VectorView<complex<float> >& b)
  {
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(b.ct() == NonConj);
    BlasTriLDivEq(A,b.Real());
    BlasTriLDivEq(A,b.Imag());
  }
#endif
#endif

  template <class T1, class T> inline void TriLDivEq(
      const GenUpperTriMatrix<T1>& A, const VectorView<T>& b)
  {
    TMVAssert(b.size() == A.size());
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
  }

  template <class T1, class T> inline void TriLDivEq(
      const GenLowerTriMatrix<T1>& A, const VectorView<T>& b)
  {
    TMVAssert(b.size() == A.size());
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
  }

  //
  // TriLDivEq M
  //

  template <class T1, class T> void RowTriLDivEq(
      const GenUpperTriMatrix<T1>& A, const MatrixView<T>& B)
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
      CVIter<T1> Aii = A.diag().begin()+N-1;
      T* Bi0 = B.ptr()+(N-1)*B.stepi();
      for(int i=N-1; i>=0; --i,--Aii,Bi0-=B.stepi()) {
	B.row(i) -= A.row(i,i+1,N) * B.Rows(i+1,N);
	if (*Aii==T1(0)) tmv_error("Singular Matrix found in UpperTriLDivEq");
	if (*Aii != T1(1)) B.row(i) /= *Aii;
      }
    }
  }

  template <class T1, class T> void ColTriLDivEq(
      const GenUpperTriMatrix<T1>& A, const MatrixView<T>& B)
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
      CVIter<T1> Ajj = A.diag().begin()+N-1;
      T* Bj0 = B.ptr()+(N-1)*B.stepi();
      for(int j=N-1; j>=0; --j,--Ajj,Bj0-=B.stepi()) {
	if (*Ajj==T1(0)) tmv_error("Singular Matrix found in UpperTriLDivEq");
	if (*Ajj != T1(1)) B.row(j) /= *Ajj;
	B.Rows(0,j) -= A.col(j,0,j) ^ B.row(j);
      }
    } 
  }

  template <class T1, class T> void RowTriLDivEq(
      const GenLowerTriMatrix<T1>& A, const MatrixView<T>& B)
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
      CVIter<T1> Aii = A.diag().begin();
      T* Bi0 = B.ptr();
      for(size_t i=0;i<N;++i,++Aii,Bi0+=B.stepi()) {
	B.row(i) -= A.row(i,0,i) * B.Rows(0,i);
	if (*Aii==T1(0)) tmv_error("Singular Matrix found in LowerTriLDivEq");
	if (*Aii != T1(1)) B.row(i) /= *Aii;
      }
    }
  }

  template <class T1, class T> void ColTriLDivEq(
      const GenLowerTriMatrix<T1>& A, const MatrixView<T>& B)
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
      CVIter<T1> Ajj = A.diag().begin();
      T* Bj0 = B.ptr();
      for(size_t j=0;j<N;++j,++Ajj,Bj0+=B.stepi()) {
	if (*Ajj==T1(0)) tmv_error("Singular Matrix found in LowerTriLDivEq");
	if (*Ajj != T1(1)) B.row(j) /= *Ajj;
	B.Rows(j+1,N) -= A.col(j,j+1,N) ^ B.row(j);
      }
    } 
  }

  template <class T1, class T> inline void NonBlasTriLDivEq(
      const GenUpperTriMatrix<T1>& A, const MatrixView<T>& B)
    // B = A^-1 * B
    // where A is a triangle matrix
  {
    //cerr<<"Upper Tri LDivEq Matrix\n";
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.colsize()>0);
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.ct() == NonConj);

    if (B.isrm()) {
      if (A.isrm()) RowTriLDivEq(A,B);
      else ColTriLDivEq(A,B);
    } else {
      for(size_t j=0;j<B.rowsize();++j) TriLDivEq(A,B.col(j));
    }
  }

  template <class T1, class T> inline void NonBlasTriLDivEq(
      const GenLowerTriMatrix<T1>& A, const MatrixView<T>& B)
    // B = A^-1 * B
    // where A is a triangle matrix
  {
    //cerr<<"Lower Tri LDivEq Matrix\n";
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.colsize()>0);
    TMVAssert(B.rowsize()>0);
    TMVAssert(B.ct() == NonConj);

    if (B.isrm()) {
      if (A.isrm()) RowTriLDivEq(A,B);
      else ColTriLDivEq(A,B);
    } else {
      for(size_t j=0;j<B.rowsize();++j) TriLDivEq(A,B.col(j));
    }
  }

#ifdef BLAS
  template <class T1, class T> inline void BlasTriLDivEq(
      const GenUpperTriMatrix<T1>& A, const MatrixView<T>& B)
  { NonBlasTriLDivEq(A,B); }
  template <class T1, class T> inline void BlasTriLDivEq(
      const GenLowerTriMatrix<T1>& A, const MatrixView<T>& B)
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
  template <class T1, class T> inline void TriLDivEq(
      const GenUpperTriMatrix<T1>& A, const MatrixView<T>& B)
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
	    UpperTriMatrix<T1,NonUnitDiag,RowMajor> tempA = A;
	    TriLDivEq(tempA,B);
	  } else {
	    UpperTriMatrix<T1,NonUnitDiag,ColMajor> tempA = A;
	    TriLDivEq(tempA,B);
	  }
	} else {
	  if (A.isrm()) {
	    UpperTriMatrix<T1,UnitDiag,RowMajor> tempA = A;
	    TriLDivEq(tempA,B);
	  } else {
	    UpperTriMatrix<T1,UnitDiag,ColMajor> tempA = A;
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
    if (Norm(BB-B0) > 0.001*Norm(B0)) {
      cerr<<"TriLDivEq: M/Upper\n";
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"Done: B = "<<B<<endl;
      cerr<<"A*B = "<<BB<<endl;
      abort();
    }
#endif
  }

  template <class T1, class T> inline void TriLDivEq(
      const GenLowerTriMatrix<T1>& A, const MatrixView<T>& B)
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
	    LowerTriMatrix<T1,NonUnitDiag,RowMajor> tempA = A;
	    TriLDivEq(tempA,B);
	  } else {
	    LowerTriMatrix<T1,NonUnitDiag,ColMajor> tempA = A;
	    TriLDivEq(tempA,B);
	  }
	} else {
	  if (A.isrm()) {
	    LowerTriMatrix<T1,UnitDiag,RowMajor> tempA = A;
	    TriLDivEq(tempA,B);
	  } else {
	    LowerTriMatrix<T1,UnitDiag,ColMajor> tempA = A;
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
    if (Norm(BB-B0) > 0.001*Norm(B0)) {
      cerr<<"TriLDivEq: M/Lower\n";
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"Done: B = "<<B<<endl;
      cerr<<"A*B = "<<BB<<endl;
      abort();
    }
#endif
  }

  template <class T1, class T> inline void RowTriLDivEq(
      const GenUpperTriMatrix<T1>& A, const UpperTriMatrixView<T>& B)
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
      CVIter<T1> Aii = A.diag().begin()+N-1;
      size_t len=1;
      size_t Bds = B.stepi()+1;
      T* Bii = B.ptr()+(N-1)*Bds;
      for(int i=N-1; i>=0; --i,--Aii,Bii-=Bds,++len) {
	B.row(i,i+1,N) -= A.row(i,i+1,N) * B.SubTriMatrix(i+1,N);
	if (*Aii==T1(0)) tmv_error("Singular Matrix found in UpperTriLDivEq");
	if (*Aii != T1(1)) B.row(i,i,N) /= *Aii;
      }
    }
  }

  template <class T1, class T> inline void ColTriLDivEq(
      const GenUpperTriMatrix<T1>& A, const UpperTriMatrixView<T>& B)
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
      CVIter<T1> Ajj = A.diag().begin()+N-1;
      size_t len=1;
      size_t Bds = B.stepi()+1;
      T* Bjj = B.ptr()+(N-1)*Bds;
      for(int j=N-1; j>=0; --j,--Ajj,Bjj-=Bds,++len) {
	if (*Ajj==T1(0)) 
	  tmv_error("Singular Matrix found in UpperTriLDivEq");
	if (*Ajj != T1(1)) B.row(j,j,N) /= *Ajj;
	B.SubMatrix(0,j,j,N) -= A.col(j,0,j) ^ B.row(j,j,N);
      }
    } 
  }

  template <class T1, class T> inline void RowTriLDivEq(
      const GenLowerTriMatrix<T1>& A, const LowerTriMatrixView<T>& B)
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
      CVIter<T1> Aii = A.diag().begin();
      T* Bi0 = B.ptr();
      for(size_t i=0;i<N;++i,++Aii,Bi0+=B.stepi()) {
	B.row(i,0,i) -= A.row(i,0,i) * B.SubTriMatrix(0,i);
	if (*Aii==T1(0)) tmv_error("Singular Matrix found in LowerTriLDivEq");
	if (*Aii != T1(1)) B.row(i,0,i+1) /= *Aii;
      }
    }
  }

  template <class T1, class T> inline void ColTriLDivEq(
      const GenLowerTriMatrix<T1>& A, const LowerTriMatrixView<T>& B)
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
      CVIter<T1> Ajj = A.diag().begin();
      T* Bj0 = B.ptr();
      for(size_t j=0;j<N;++j,++Ajj,Bj0+=B.stepi()) {
	if (*Ajj==T1(0)) 
	  tmv_error("Singular Matrix found in LowerTriLDivEq");
	if (*Ajj != T1(1)) B.row(j,0,j+1) /= *Ajj;
	B.SubMatrix(j+1,N,0,j+1) -= A.col(j,j+1,N) ^ B.row(j,0,j+1);
      }
    }
  }

  template <class T1, class T> inline void TriLDivEq(
      const GenUpperTriMatrix<T1>& A, const UpperTriMatrixView<T>& B)
    // B = A^-1 * B
    // where A is a triangle matrix
  {
#ifdef XDEBUG
    UpperTriMatrix<T> B0 = B;
#endif
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || A.isunit());

    if (B.isconj()) TriLDivEq(A.Conjugate(),B.Conjugate());
    else if (B.size() > 0) 
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

#ifdef XDEBUG
    UpperTriMatrix<T> BB = A*B;
    if (Norm(BB-B0) > 0.001*Norm(B0)) {
      cerr<<"TriLDivEq: Upper/Upper\n";
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"Done: B = "<<B<<endl;
      cerr<<"A*B = "<<BB<<endl;
      abort();
    }
#endif
  }

  template <class T1, class T> inline void TriLDivEq(
      const GenLowerTriMatrix<T1>& A, const LowerTriMatrixView<T>& B)
    // B = A^-1 * B
    // where A is a triangle matrix
  {
#ifdef XDEBUG
    LowerTriMatrix<T> B0 = B;
#endif

    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || A.isunit());

    if (B.isconj()) TriLDivEq(A.Conjugate(),B.Conjugate());
    else if (B.size() > 0) 
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

#ifdef XDEBUG
    LowerTriMatrix<T> BB = A*B;
    if (Norm(BB-B0) > 0.001*Norm(B0)) {
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


