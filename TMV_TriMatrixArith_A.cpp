
#include "TMV_VectorArith_Inline.h"
#include "TMV_Tri.h"

namespace tmv {

  //
  // MultMV
  //

  template <class T, class Ta> void RowMajorMultEqMV(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& x)
  {
    TMVAssert(A.isrm());
    TMVAssert(x.step()==1);
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);

    size_t N = x.size();
    while (N>0 && (x(N-1) == T(0))) --N;
    if (N == 0) return;

    const int diagstep = A.stepi()+1;
    VIt<T,Unit,NonConj> xi = x.begin();
    const Ta* Aptr = A.cptr();
    if (A.isunit()) {
      Aptr++;
      for(size_t i=0,len=N-1; len>0; ++i,--len,++xi,Aptr+=diagstep) {
	T temp;
	if (A.isconj())
	  DoMultVV(temp,CVIt<Ta,Unit,Conj>(Aptr,1),
	      CVIt<T,Unit,NonConj>(xi+1),len);
	else
	  DoMultVV(temp,CVIt<Ta,Unit,NonConj>(Aptr,1),
	      CVIt<T,Unit,NonConj>(xi+1),len);
	*xi += temp;
      }
    }
    else {
      for(size_t i=0,len=N; len>0; ++i,--len,++xi,Aptr+=diagstep) {
	T temp;
	if (A.isconj())
	  DoMultVV(temp,CVIt<Ta,Unit,Conj>(Aptr,1),
	      CVIt<T,Unit,NonConj>(xi),len);
	else
	  DoMultVV(temp,CVIt<Ta,Unit,NonConj>(Aptr,1),
	      CVIt<T,Unit,NonConj>(xi),len);
	*xi = temp;
      }
    }
  }

  template <class T, class Ta> void RowMultEqMV(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& x)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);

    size_t N = x.size();
    while (N>0 && (x(N-1) == T(0))) --N;
    if (N == 0) return;

    VIt<T,Step,NonConj> xit = x.begin();
    if (A.isunit()) {
      for(size_t i=0; i<N; ++i,++xit) 
	*xit += A.row(i,i+1,N) * x.SubVector(i+1,N);
    }
    else {
      for(size_t i=0; i<N; ++i,++xit) 
	*xit = A.row(i,i,N) * x.SubVector(i,N);
    }
  }

  template <class T, class Ta> void ColMajorMultEqMV(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& x)
  {
    TMVAssert(A.iscm());
    TMVAssert(x.step()==1);
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);

    size_t N = x.size();
    while (N>0 && (x(N-1) == T(0))) --N;
    if (N == 0) return;

    VIt<T,Unit,NonConj> x0 = x.begin();
    VIt<T,Unit,NonConj> xj = x.begin()+1;
    const Ta* Aptr = A.cptr()+A.stepj();
    if (A.isunit()) {
      for(size_t j=1; j<N; ++j,++xj,Aptr+=A.stepj()) {
	if (*xj != T(0)) {
	  if (A.isconj())
	    DoAddVV(*xj,CVIt<Ta,Unit,Conj>(Aptr,1),x0,j);
	  else
	    DoAddVV(*xj,CVIt<Ta,Unit,NonConj>(Aptr,1),x0,j);
	}
      }
    }
    else {
      CVIter<Ta> Ajj = A.diag().begin();
      *x0 *= *Ajj;  ++Ajj;
      for(size_t j=1; j<N; ++j,++xj,++Ajj,Aptr+=A.stepj()) {
	if (*xj != T(0)) {
	  if (A.isconj())
	    DoAddVV(*xj,CVIt<Ta,Unit,Conj>(Aptr,1),x0,j);
	  else
	    DoAddVV(*xj,CVIt<Ta,Unit,NonConj>(Aptr,1),x0,j);
	}
	*xj *= *Ajj;
      }
    }
  }

  template <class T, class Ta> void RowMajorMultEqMV(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& x)
  {
    TMVAssert(A.isrm());
    TMVAssert(x.step()==1);
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);

    const int N = x.size();
    int i1 = 0;
    while (i1<N && (x(i1) == T(0))) ++i1;
    if (i1 == N) return;

    CVIt<T,Unit,NonConj> xi1 = x.begin()+i1;
    VIt<T,Unit,NonConj> xi = x.begin()+N-1;
    const Ta* Aptr = A.cptr()+(N-1)*A.stepi()+i1;
    if (A.isunit()) {
      for(int i=N-1,len=i-i1; len>0; --i,--len,--xi,Aptr-=A.stepi()) {
	T temp;
	if (A.isconj())
	  DoMultVV(temp,CVIt<Ta,Unit,Conj>(Aptr,1),xi1,len);
	else
	  DoMultVV(temp,CVIt<Ta,Unit,NonConj>(Aptr,1),xi1,len);
	*xi += temp;
      }
    }
    else {
      for(int i=N-1,len=i+1-i1; len>0; --i,--len,--xi,Aptr-=A.stepi()) {
	T temp;
	if (A.isconj())
	  DoMultVV(temp,CVIt<Ta,Unit,Conj>(Aptr,1),xi1,len);
	else
	  DoMultVV(temp,CVIt<Ta,Unit,NonConj>(Aptr,1),xi1,len);
	*xi = temp;
      }
    }
  }

  template <class T, class Ta> void RowMultEqMV(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& x)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);

    const int N = x.size();
    int i1 = 0;
    while (i1<N && (x(i1) == T(0))) ++i1;
    if (i1 == N) return;

    VIt<T,Step,NonConj> xit = x.begin()+N-1;
    if (A.isunit()) {
      for(int i=N-1; i>=i1; --i,--xit) 
	*xit += A.row(i,i1,i) * x.SubVector(i1,i);
    }
    else {
      for(int i=N-1; i>=i1; --i,--xit) 
	*xit = A.row(i,i1,i+1) * x.SubVector(i1,i+1);
    }
  }

  template <class T, class Ta> void ColMajorMultEqMV(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& x)
  {
    TMVAssert(A.iscm());
    TMVAssert(x.step()==1);
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);

    const int N = x.size();
    int i1 = 0;
    while (i1<N && (x(i1) == T(0))) ++i1;
    if (i1 == N) return;

    VIt<T,Unit,NonConj> xj = x.begin()+N-2;
    const int diagstep = A.stepj()+1;
    const Ta* Aptr = A.cptr()+(N-2)*diagstep+1;
    if (A.isunit()) {
      for(int j=N-2,len=1; j>=i1; --j,++len,--xj,Aptr-=diagstep) {
	if (*xj != T(0)) {
	  if (A.isconj())
	    DoAddVV(*xj,CVIt<Ta,Unit,Conj>(Aptr,1),xj+1,len);
	  else
	    DoAddVV(*xj,CVIt<Ta,Unit,NonConj>(Aptr,1),xj+1,len);
	}
      }
    }
    else {
      CVIter<Ta> Ajj = A.diag().begin()+N-1;
      *(xj+1) *= *Ajj;  --Ajj;
      for(int j=N-2,len=1; j>=i1; --j,++len,--xj,--Ajj,Aptr-=diagstep) {
	if (*xj != T(0)) {
	  if (A.isconj())
	    DoAddVV(*xj,CVIt<Ta,Unit,Conj>(Aptr,1),xj+1,len);
	  else
	    DoAddVV(*xj,CVIt<Ta,Unit,NonConj>(Aptr,1),xj+1,len);
	}
	*xj *= *Ajj;
      }
    }
  }

  template <class T, class Ta> inline void NonBlasMultEqMV(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& x)
    // x = A * x
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);

    if (x.step()==1) {
      if (A.isrm()) RowMajorMultEqMV(A,x);
      else if (A.iscm()) ColMajorMultEqMV(A,x);
      else RowMultEqMV(A,x);
    }
    else RowMultEqMV(A,x);
  }

  template <class T, class Ta> inline void NonBlasMultEqMV(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& x)
    // x = A * x
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);

    if (x.step()==1) {
      if (A.isrm()) RowMajorMultEqMV(A,x);
      else if (A.iscm()) ColMajorMultEqMV(A,x);
      else RowMultEqMV(A,x);
    }
    else RowMultEqMV(A,x);
  }

#ifdef BLAS
  template <class T, class Ta> inline void BlasMultEqMV(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& x)
  { NonBlasMultEqMV(A,x); }
  template <class T, class Ta> inline void BlasMultEqMV(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& x)
  { NonBlasMultEqMV(A,x); }
  inline void BlasMultEqMV( 
      const GenUpperTriMatrix<double>& A, const VectorView<double>& x)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    cblas_dtrmv(A.isrm() ? CblasRowMajor : CblasColMajor,
	CblasUpper, CblasNoTrans, A.isunit() ? CblasUnit : CblasNonUnit,
	A.size(),A.cptr(), A.isrm() ? A.stepi() : A.stepj(),
	x.ptr(),x.step()); 
  }
  inline void BlasMultEqMV( 
      const GenLowerTriMatrix<double>& A, const VectorView<double>& x)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    cblas_dtrmv(A.isrm() ? CblasRowMajor : CblasColMajor,
	CblasLower, CblasNoTrans, A.isunit() ? CblasUnit : CblasNonUnit,
	A.size(),A.cptr(), A.isrm() ? A.stepi() : A.stepj(),
	x.ptr(),x.step()); 
  }
  inline void BlasMultEqMV(const GenUpperTriMatrix<complex<double> >& A,
      const VectorView<complex<double> >& x)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);
    CBLAS_ORDER order = (A.iscm() == A.isconj()) ? 
      CblasRowMajor : CblasColMajor;
    CBLAS_UPLO uplo = A.isconj() ? CblasLower : CblasUpper;
    CBLAS_TRANSPOSE tran = A.isconj() ?  CblasConjTrans : CblasNoTrans;
    cblas_ztrmv(order, uplo, tran, A.isunit() ? CblasUnit : CblasNonUnit,
	A.size(),A.cptr(), A.isrm() ? A.stepi() : A.stepj(),
	x.ptr(),x.step()); 
  }
  inline void BlasMultEqMV(const GenLowerTriMatrix<complex<double> >& A,
      const VectorView<complex<double> >& x)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);
    CBLAS_ORDER order = (A.iscm() == A.isconj()) ? 
      CblasRowMajor : CblasColMajor;
    CBLAS_UPLO uplo = A.isconj() ? CblasUpper : CblasLower;
    CBLAS_TRANSPOSE tran = A.isconj() ?  CblasConjTrans : CblasNoTrans;
    cblas_ztrmv(order, uplo, tran, A.isunit() ? CblasUnit : CblasNonUnit,
	A.size(),A.cptr(), A.isrm() ? A.stepi() : A.stepj(),
	x.ptr(),x.step()); 
  }
  inline void BlasMultEqMV(const GenUpperTriMatrix<double>& A,
      const VectorView<complex<double> >& x)
  {
    TMVAssert(x.ct() == NonConj);
    BlasMultEqMV(A,x.Real());
    BlasMultEqMV(A,x.Imag());
  }
  inline void BlasMultEqMV(const GenLowerTriMatrix<double>& A,
      const VectorView<complex<double> >& x)
  {
    TMVAssert(x.ct() == NonConj);
    BlasMultEqMV(A,x.Real());
    BlasMultEqMV(A,x.Imag());
  }
#ifndef NOFLOAT
  inline void BlasMultEqMV( 
      const GenUpperTriMatrix<float>& A, const VectorView<float>& x)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    cblas_strmv(A.isrm() ? CblasRowMajor : CblasColMajor,
	CblasUpper, CblasNoTrans, A.isunit() ? CblasUnit : CblasNonUnit,
	A.size(),A.cptr(), A.isrm() ? A.stepi() : A.stepj(),
	x.ptr(),x.step()); 
  }
  inline void BlasMultEqMV( 
      const GenLowerTriMatrix<float>& A, const VectorView<float>& x)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    cblas_strmv(A.isrm() ? CblasRowMajor : CblasColMajor,
	CblasLower, CblasNoTrans, A.isunit() ? CblasUnit : CblasNonUnit,
	A.size(),A.cptr(), A.isrm() ? A.stepi() : A.stepj(),
	x.ptr(),x.step()); 
  }
  inline void BlasMultEqMV(const GenUpperTriMatrix<complex<float> >& A,
      const VectorView<complex<float> >& x)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);
    CBLAS_ORDER order = (A.iscm() == A.isconj()) ? 
      CblasRowMajor : CblasColMajor;
    CBLAS_UPLO uplo = A.isconj() ? CblasLower : CblasUpper;
    CBLAS_TRANSPOSE tran = A.isconj() ? CblasConjTrans : CblasNoTrans;
    cblas_ctrmv(order, uplo, tran, A.isunit() ? CblasUnit : CblasNonUnit,
	A.size(),A.cptr(), A.isrm() ? A.stepi() : A.stepj(),
	x.ptr(),x.step()); 
  }
  inline void BlasMultEqMV(const GenLowerTriMatrix<complex<float> >& A,
      const VectorView<complex<float> >& x)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);
    CBLAS_ORDER order = (A.iscm() == A.isconj()) ? 
      CblasRowMajor : CblasColMajor;
    CBLAS_UPLO uplo = A.isconj() ? CblasUpper : CblasLower;
    CBLAS_TRANSPOSE tran = A.isconj() ?  CblasConjTrans : CblasNoTrans;
    cblas_ctrmv(order, uplo, tran, A.isunit() ? CblasUnit : CblasNonUnit,
	A.size(),A.cptr(), A.isrm() ? A.stepi() : A.stepj(),
	x.ptr(),x.step()); 
  }
  inline void BlasMultEqMV(const GenUpperTriMatrix<float>& A,
      const VectorView<complex<float> >& x)
  {
    TMVAssert(x.ct() == NonConj);
    BlasMultEqMV(A,x.Real());
    BlasMultEqMV(A,x.Imag());
  }
  inline void BlasMultEqMV(const GenLowerTriMatrix<float>& A,
      const VectorView<complex<float> >& x)
  {
    TMVAssert(x.ct() == NonConj);
    BlasMultEqMV(A,x.Real());
    BlasMultEqMV(A,x.Imag());
  }
#endif // NOFLOAT
#endif // BLAS

  template <class T, class Ta> inline void MultEqMV(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& x)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    if (x.isconj()) MultEqMV(A.Conjugate(),x.Conjugate());
    else 
#ifdef BLAS
      if ( (A.isrm() || A.iscm()) && x.step()>0)
	BlasMultEqMV(A,x);
      else
#endif
	NonBlasMultEqMV(A,x);
  }

  template <class T, class Ta> inline void MultEqMV(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& x)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    if (x.isconj()) MultEqMV(A.Conjugate(),x.Conjugate());
    else 
#ifdef BLAS
      if ( (A.isrm() || A.iscm()) && x.step()>0)
	BlasMultEqMV(A,x);
      else
#endif
	NonBlasMultEqMV(A,x);
  }

  template <class T, class Ta, class Tx>
    void MultMV(const T alpha, const GenUpperTriMatrix<Ta>& A,
	const GenVector<Tx>& x, const T beta, const VectorView<T>& y)
    // y = alpha * A * x + beta * y    
    { 
      TMVAssert(A.size() == x.size());
      TMVAssert(A.size() == y.size());

      if (y.size() > 0) {
	if (alpha==T(0)) {
	  if (beta != T(1)) y *= beta;
	} else {
	  if (beta == T(0)) {
	    y = x;
	    MultEqMV(A,y);
	    if(alpha != T(1)) y *= alpha;
	  } else {
	    Vector<T> temp = x;
	    MultEqMV(A,temp.View());
	    if (beta != T(1)) y *= beta;
	    y += alpha * temp;
	  }
	} 
      }
    }

  template <class T, class Ta, class Tx>
    void MultMV(const T alpha, const GenLowerTriMatrix<Ta>& A,
	const GenVector<Tx>& x, const T beta, const VectorView<T>& y)
    // y = alpha * A * x + beta * y    
    { 
      TMVAssert(A.size() == x.size());
      TMVAssert(A.size() == y.size());

      if (y.size() > 0) {
	if (alpha==T(0)) {
	  if (beta != T(1)) y *= beta;
	} else {
	  if (beta == T(0)) {
	    y = x;
	    MultEqMV(A,y);
	    if(alpha != T(1)) y *= alpha;
	  } else {
	    Vector<T> temp = x;
	    MultEqMV(A,temp.View());
	    if (beta != T(1)) y *= beta;
	    y += alpha * temp;
	  }
	} 
      }
    }

#define InstFile "TMV_TriMatrixArith_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


