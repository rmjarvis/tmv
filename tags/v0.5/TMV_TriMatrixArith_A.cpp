
#include "TMV_VectorArith_Inline.h"
#include "TMV_Tri.h"

//#define XDEBUG

namespace tmv {

  //
  // MultMV
  //

  template <class T, class Ta> void RowMajorMultEqMV(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& x)
  {
    //cerr<<"RowMajorMultEqMV Upper\n";
    TMVAssert(A.isrm());
    TMVAssert(x.step()==1);
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);

    size_t N = x.size();
    while (N>0 && (x(N-1) == T(0))) --N;
    if (N == 0) return;

    const int ds = A.stepi()+1;
    VIt<T,Unit,NonConj> xi = x.begin();
    const Ta* Aptr = A.cptr();
    size_t len = N;
    if (A.isunit()) { ++Aptr, --len; }
    for(size_t i=0; len>0; ++i,--len,++xi,Aptr+=ds) {
      // temp = A.row(i,ii,N) * x.SubVector(ii,N);
      T temp;
      CVIt<T,Unit,NonConj> xii = A.isunit() ? xi+1 : xi;
      if (A.isconj())
	DoMultVV(temp,CVIt<Ta,Unit,Conj>(Aptr,1),xii,len);
      else
	DoMultVV(temp,CVIt<Ta,Unit,NonConj>(Aptr,1),xii,len);
      if (A.isunit()) *xi += temp;
      else *xi = temp;
    }
  }

  template <class T, class Ta> void RowMultEqMV(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& x)
  {
    //cerr<<"RowMultEqMV Upper\n";
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);

    size_t N = x.size();
    while (N>0 && (x(N-1) == T(0))) --N;
    if (N == 0) return;

    VIt<T,Step,NonConj> xit = x.begin();
    for(size_t i=0; i<N; ++i,++xit) {
      if (A.isunit()) *xit += A.row(i,i+1,N) * x.SubVector(i+1,N);
      else *xit = A.row(i,i,N) * x.SubVector(i,N);
    }
  }

  template <class T, class Ta> void ColMajorMultEqMV(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& x)
  {
    //cerr<<"ColMajorMultEqMV Upper\n";
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
	// x.SubVector(0,j) += *xj * A.col(j,0,j);
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
	// x.SubVector(0,j) += *xj * A.col(j,0,j);
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
    //cerr<<"RowMajorMultEqMV Lower\n";
    TMVAssert(A.isrm());
    TMVAssert(x.step()==1);
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);

    const size_t N = x.size();
    size_t i1 = 0;
    while (i1<N && (x(i1) == T(0))) ++i1;
    if (i1 == N) return;

    CVIt<T,Unit,NonConj> xi1 = x.begin()+i1;
    VIt<T,Unit,NonConj> xi = x.begin()+N-1;
    const Ta* Aptr = A.cptr()+(N-1)*A.stepi()+i1;
    size_t len = N-i1;
    if (A.isunit()) --len;
    for(size_t i=N-1; len>0; --i,--len,--xi,Aptr-=A.stepi()) {
      // temp = A.row(i,i1,ii) * x.SubVector(i1,ii);
      T temp;
      if (A.isconj())
	DoMultVV(temp,CVIt<Ta,Unit,Conj>(Aptr,1),xi1,len);
      else
	DoMultVV(temp,CVIt<Ta,Unit,NonConj>(Aptr,1),xi1,len);
      if (A.isunit()) *xi += temp;
      else *xi = temp;
    }
  }

  template <class T, class Ta> void RowMultEqMV(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& x)
  {
    //cerr<<"RowMultEqMV Lower\n";
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);

    const int N = x.size();
    int i1 = 0;
    while (i1<N && (x(i1) == T(0))) ++i1;
    if (i1 == N) return;

    VIt<T,Step,NonConj> xit = x.begin()+N-1;
    for(int i=N-1; i>=i1; --i,--xit) {
      if (A.isunit()) 
	*xit += A.row(i,i1,i) * x.SubVector(i1,i);
      else 
	*xit = A.row(i,i1,i+1) * x.SubVector(i1,i+1);
    }
  }

  template <class T, class Ta> void ColMajorMultEqMV(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& x)
  {
    //cerr<<"ColMajorMultEqMV Lower\n";
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
    const int ds = A.stepj()+1;
    const Ta* Aptr = A.cptr()+(N-2)*ds+1;
    if (A.isunit()) {
      for(int j=N-2,len=1; j>=i1; --j,++len,--xj,Aptr-=ds) {
	// x.SubVector(j+1,N) += *xj * A.col(j,j+1,N);
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
      for(int j=N-2,len=1; j>=i1; --j,++len,--xj,--Ajj,Aptr-=ds) {
	// x.SubVector(j+1,N) += *xj * A.col(j,j+1,N);
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
#ifdef XDEBUG
    Vector<T> x2 = Matrix<T>(A) * x;
    Vector<T> x0 = x;
#endif
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
#ifdef XDEBUG
    if (Norm(x-x2) > 0.001*Norm(x)) {
      cerr<<"MultEqMV: \n";
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"x = "<<Type(x)<<" step "<<x.step()<<"  "<<x0<<endl;
      cerr<<"-> x = "<<x<<endl;
      cerr<<"x2 = "<<x2<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta> inline void MultEqMV(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& x)
  {
#ifdef XDEBUG
    Vector<T> x2 = Matrix<T>(A) * x;
    Vector<T> x0 = x;
#endif
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
#ifdef XDEBUG
    if (Norm(x-x2) > 0.001*Norm(x)) {
      cerr<<"MultEqMV: \n";
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"x = "<<Type(x)<<" step "<<x.step()<<"  "<<x0<<endl;
      cerr<<"-> x = "<<x<<endl;
      cerr<<"x2 = "<<x2<<endl;
      abort();
    }
#endif
  }

  template <class T1, class Ta, class Tx, class Ty> void RowMajorAddMultMV(
      const T1 alpha, const GenUpperTriMatrix<Ta>& A,
      const GenVector<Tx>& x, const VectorView<Ty>& y)
  {
    //cerr<<"RowMajor AddMult Upper\n";
    TMVAssert(A.isrm());
    TMVAssert(x.step()==1);
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(alpha != T1(0));
    TMVAssert(y.ct() == NonConj);

    size_t N = x.size();
    while (N>0 && (x(N-1) == Tx(0))) --N;
    if (N == 0) return;

    const int ds = A.stepi()+1;
    VIt<Ty,Step,NonConj> yi = y.begin();
    CVIter<Tx> xi = x.begin();
    const Ta* Aptr = A.cptr();
    size_t len = N;
    if (A.isunit()) { ++Aptr; --len; }

    for(size_t i=0; len>0; ++i,--len,++xi,++yi,Aptr+=ds) {
      Ty temp;
      if (x.isconj()) {
	CVIt<Tx,Unit,Conj> xii = A.isunit() ? xi+1 : xi;
	if (A.isconj())
	  DoMultVV(temp,CVIt<Ta,Unit,Conj>(Aptr,1),xii,len);
	else
	  DoMultVV(temp,CVIt<Ta,Unit,NonConj>(Aptr,1),xii,len);
      } else {
	CVIt<Tx,Unit,NonConj> xii = A.isunit() ? xi+1 : xi;
	if (A.isconj()) 
	  DoMultVV(temp,CVIt<Ta,Unit,Conj>(Aptr,1),xii,len);
	else
	  DoMultVV(temp,CVIt<Ta,Unit,NonConj>(Aptr,1),xii,len);
      }
      if (A.isunit()) temp += *xi;
      if (alpha != T1(1)) {
	if (alpha == T1(-1)) temp = -temp;
	else temp *= alpha;
      }
      *yi += temp;
    }
    if (A.isunit()) {
      if (alpha == T1(1)) *yi += *xi;
      else if (alpha == T1(-1)) *yi -= *xi;
      else *yi += *xi * alpha;
    }
  }

  template <class T, class Ta, class Tx> void RowAddMultMV(
      const T alpha, const GenUpperTriMatrix<Ta>& A,
      const GenVector<Tx>& x, const VectorView<T>& y)
  {
    //cerr<<"Row AddMult Upper\n";
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(y.ct() == NonConj);

    size_t N = x.size();
    while (N>0 && (x(N-1) == T(0))) --N;
    if (N == 0) return;

    VIt<T,Step,NonConj> yit = y.begin();
    if (A.isunit()) {
      CVIter<Tx> xit = x.begin();
      for(size_t i=0; i<N; ++i,++xit,++yit) 
	if (alpha == T(1))
	  *yit += *xit + A.row(i,i+1,N) * x.SubVector(i+1,N);
	else
	  *yit += alpha * (*xit + A.row(i,i+1,N) * x.SubVector(i+1,N));
    }
    else {
      for(size_t i=0; i<N; ++i,++yit) 
	if (alpha == T(1)) *yit += A.row(i,i,N) * x.SubVector(i,N);
	else *yit += alpha * A.row(i,i,N) * x.SubVector(i,N);
    }
  }

  template <class T1, class Ta, class Tx, class Ty> void ColMajorAddMultMV(
      const T1 alpha, const GenUpperTriMatrix<Ta>& A,
      const GenVector<Tx>& x, const VectorView<Ty>& y)
  {
    //cerr<<"ColMajor AddMult Upper\n";
    TMVAssert(A.iscm());
    TMVAssert(y.step()==1);
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != T1(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.ct() == NonConj);

    size_t N = x.size();
    while (N>0 && (x(N-1) == Tx(0))) --N;
    if (N == 0) return;

    VIt<Ty,Unit,NonConj> y0 = y.begin();
    CVIter<Tx> xj = x.begin();
    const Ta* Aptr = A.cptr();
    if (A.isunit()) {
      VIt<Ty,Unit,NonConj> yj = y0;
      if (*xj != Tx(0)) {
	if (alpha == T1(1)) *yj += *xj*alpha;
	else *yj += *xj*alpha;
      }
      ++xj; ++yj; Aptr+=A.stepj();
      for(size_t j=1; j<N; ++j,++xj,++yj,Aptr+=A.stepj()) {
	if (*xj != Tx(0)) {
	  Ty ax = *xj;
	  if (alpha != T1(1)) ax *= alpha;
	  // y.SubVector(0,j) += ax * A.col(j,0,j);
	  if (IMAG(ax) == RealType(Ty)(0)) {
	    if (A.isconj())
	      DoAddVV(REAL(ax),CVIt<Ta,Unit,Conj>(Aptr,1),y0,j);
	    else
	      DoAddVV(REAL(ax),CVIt<Ta,Unit,NonConj>(Aptr,1),y0,j);
	    *yj += REAL(ax);
	  } else {
	    if (A.isconj())
	      DoAddVV(ax,CVIt<Ta,Unit,Conj>(Aptr,1),y0,j);
	    else
	      DoAddVV(ax,CVIt<Ta,Unit,NonConj>(Aptr,1),y0,j);
	    *yj += ax;
	  }
	}
      }
    }
    else {
      for(size_t j=0; j<N; ++j,++xj,Aptr+=A.stepj()) {
	if (*xj != Tx(0)) {
	  Ty ax = *xj;
	  if (alpha != T1(1)) ax *= alpha;
	  // y.SubVector(0,j+1) += ax * A.col(j,0,j+1);
	  if (IMAG(ax) == RealType(Ty)(0))
	    if (A.isconj())
	      DoAddVV(REAL(ax),CVIt<Ta,Unit,Conj>(Aptr,1),y0,j+1);
	    else
	      DoAddVV(REAL(ax),CVIt<Ta,Unit,NonConj>(Aptr,1),y0,j+1);
	  else
	    if (A.isconj())
	      DoAddVV(ax,CVIt<Ta,Unit,Conj>(Aptr,1),y0,j+1);
	    else
	      DoAddVV(ax,CVIt<Ta,Unit,NonConj>(Aptr,1),y0,j+1);
	}
      }
    }
  }

  template <class T1, class Ta, class Tx, class Ty> void RowMajorAddMultMV(
      const T1 alpha, const GenLowerTriMatrix<Ta>& A,
      const GenVector<Tx>& x, const VectorView<Ty>& y)
  {
    //cerr<<"RowMajor AddMult Lower\n";
    TMVAssert(A.isrm());
    TMVAssert(x.step()==1);
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(alpha != T1(0));
    TMVAssert(y.ct() == NonConj);

    const size_t N = x.size();
    size_t i1 = 0;
    while (i1<N && (x(i1) == Tx(0))) ++i1;
    if (i1 == N) return;

    const Tx* xptr = x.cptr()+i1;
    VIt<Ty,Step,NonConj> yi = y.begin()+i1;
    const Ta* Aptr = A.cptr()+i1*(A.stepi()+1);
    if (A.isunit()) {
      CVIter<Tx> xi = x.begin()+i1;
      if (alpha == T1(1)) *yi += *xi;
      else if (alpha == T1(-1)) *yi -= *xi;
      else *yi += alpha * (*xi);
      ++xi,++yi,Aptr+=A.stepi();
      for(size_t i=i1+1,len=1; i<N; ++i,++xi,++yi,++len,Aptr+=A.stepi()) {
	// temp = A.row(i,i1,i) * x.SubVector(i1,i);
	Ty temp;
	if (A.isconj())
	  if (x.isconj())
	    DoMultVV(temp,CVIt<Ta,Unit,Conj>(Aptr,1),
		CVIt<Tx,Unit,Conj>(xptr,1),len);
	  else
	    DoMultVV(temp,CVIt<Ta,Unit,Conj>(Aptr,1),
		CVIt<Tx,Unit,NonConj>(xptr,1),len);
	else
	  if (x.isconj())
	    DoMultVV(temp,CVIt<Ta,Unit,NonConj>(Aptr,1),
		CVIt<Tx,Unit,Conj>(xptr,1),len);
	  else
	    DoMultVV(temp,CVIt<Ta,Unit,NonConj>(Aptr,1),
		CVIt<Tx,Unit,NonConj>(xptr,1),len);
	if (alpha == T1(1)) *yi += *xi + temp;
	else if (alpha == T1(-1)) *yi -= *xi + temp;
	else *yi += alpha * (*xi + temp);
      }
    }
    else {
      for(size_t i=i1,len=1; i<N; ++i,++len,++yi,Aptr+=A.stepi()) {
	// temp = A.row(i,i1,i+1) * x.SubVector(i1,i+1);
	Ty temp;
	if (A.isconj())
	  if (x.isconj())
	    DoMultVV(temp,CVIt<Ta,Unit,Conj>(Aptr,1),
		CVIt<Tx,Unit,Conj>(xptr,1),len);
	  else
	    DoMultVV(temp,CVIt<Ta,Unit,Conj>(Aptr,1),
		CVIt<Tx,Unit,NonConj>(xptr,1),len);
	else
	  if (x.isconj())
	    DoMultVV(temp,CVIt<Ta,Unit,NonConj>(Aptr,1),
		CVIt<Tx,Unit,Conj>(xptr,1),len);
	  else
	    DoMultVV(temp,CVIt<Ta,Unit,NonConj>(Aptr,1),
		CVIt<Tx,Unit,NonConj>(xptr,1),len);
	if (alpha == T1(1)) *yi += temp;
	else if (alpha == T1(-1)) *yi -= temp;
	else *yi += alpha * temp;
      }
    }
  }

  template <class T, class Ta, class Tx> void RowAddMultMV(
      const T alpha, const GenLowerTriMatrix<Ta>& A,
      const GenVector<Tx>& x, const VectorView<T>& y)
  {
    //cerr<<"Row AddMult Lower\n";
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.ct() == NonConj);

    const size_t N = x.size();
    size_t i1 = 0;
    while (i1<N && (x(i1) == Tx(0))) ++i1;
    if (i1 == N) return;

    VIt<T,Step,NonConj> yit = y.begin()+i1;
    if (A.isunit()) {
      CVIter<Tx> xit = x.begin()+i1;
      for(size_t i=i1; i<N; ++i,++xit,++yit) 
	if (alpha == T(1)) 
	  *yit += *xit + A.row(i,i1,i) * x.SubVector(i1,i);
	else 
	  *yit += alpha*(*xit + A.row(i,i1,i) * x.SubVector(i1,i));
    }
    else {
      for(size_t i=i1; i<N; ++i,++yit) 
	if (alpha == T(1)) 
	  *yit += A.row(i,i1,i+1) * x.SubVector(i1,i+1);
	else 
	  *yit += alpha * A.row(i,i1,i+1) * x.SubVector(i1,i+1);
    }
  }

  template <class T1, class Ta, class Tx, class Ty> void ColMajorAddMultMV(
      const T1 alpha, const GenLowerTriMatrix<Ta>& A,
      const GenVector<Tx>& x, const VectorView<Ty>& y)
  {
    //cerr<<"ColMajor AddMult Lower\n";
    TMVAssert(A.iscm());
    TMVAssert(y.step()==1);
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != T1(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.ct() == NonConj);

    const size_t N = x.size();
    size_t i1 = 0;
    while (i1<N && (x(i1) == Tx(0))) ++i1;
    if (i1 == N) return;

    const int ds = A.stepj()+1;
    VIt<Ty,Unit,NonConj> yj = y.begin()+i1;
    CVIter<Tx> xj = x.begin()+i1;
    const Ta* Aptr = A.cptr()+i1*ds;
    size_t len = N-i1;
    if (A.isunit()) { ++Aptr; --len; }

    for(size_t j=i1; len>0; ++j,--len,++xj,++yj,Aptr+=ds) {
      if (*xj != Tx(0)) {
	Ty ax = *xj;
	if (alpha != T1(1)) ax *= alpha;
	VIt<Ty,Unit,NonConj> yjj = A.isunit() ? yj+1 : yj;
	// y.SubVector(jj,N) += ax * A.col(j,jj,N);
	if (IMAG(ax) == RealType(Ty)(0)) {
	  if (A.isconj())
	    DoAddVV(REAL(ax),CVIt<Ta,Unit,Conj>(Aptr,1),yjj,len);
	  else
	    DoAddVV(REAL(ax),CVIt<Ta,Unit,NonConj>(Aptr,1),yjj,len);
	  if (A.isunit()) *yj += REAL(ax);
	}
	else {
	  if (A.isconj())
	    DoAddVV(ax,CVIt<Ta,Unit,Conj>(Aptr,1),yjj,len);
	  else
	    DoAddVV(ax,CVIt<Ta,Unit,NonConj>(Aptr,1),yjj,len);
	  if (A.isunit()) *yj += ax;
	}
      }
    }
    if (A.isunit() && *xj != Tx(0)) {
      if (alpha == T1(1)) *yj += *xj;
      else *yj += *xj*alpha;
    }
  }

  template <class T, class Ta, class Tx> void AddMultMV(
      const T alpha, const GenUpperTriMatrix<Ta>& A,
      const GenVector<Tx>& x, const VectorView<T>& y)
    // y += alpha * A * x
  {
#ifdef XDEBUG
    Vector<T> y2 = y + alpha *Matrix<T>(A)*x;
    Vector<T> y0 = y;
#endif
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());

    if (y.size() > 0) {
      if (y.isconj()) AddMultMV(CONJ(alpha),A.Conjugate(),
	  x.Conjugate(),y.Conjugate());
      else if (A.isrm() && x.step()==1) 
	if (IMAG(alpha) == RealType(T)(0)) 
	  RowMajorAddMultMV(REAL(alpha),A,x,y);
	else RowMajorAddMultMV(alpha,A,x,y);
      else if (A.iscm() && y.step()==1) 
	if (IMAG(alpha) == RealType(T)(0)) 
	  ColMajorAddMultMV(REAL(alpha),A,x,y);
	else ColMajorAddMultMV(alpha,A,x,y);
      else RowAddMultMV(alpha,A,x,y);
    }
#ifdef XDEBUG
    if (Norm(y-y2) > 0.001*Norm(x)) {
      cerr<<"AddMultMV: alpha = "<<alpha<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"x = "<<Type(x)<<" step "<<x.step()<<"  "<<x<<endl;
      cerr<<"y = "<<Type(y)<<" step "<<y.step()<<"  "<<y0<<endl;
      cerr<<"-> y = "<<y<<endl;
      cerr<<"y2 = "<<y2<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta, class Tx> void AddMultMV(
      const T alpha, const GenLowerTriMatrix<Ta>& A,
      const GenVector<Tx>& x, const VectorView<T>& y)
    // y += alpha * A * x
  {
#ifdef XDEBUG
    Vector<T> y2 = y + alpha *Matrix<T>(A)*x;
    Vector<T> y0 = y;
#endif
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());

    if (y.size() > 0) {
      if (y.isconj()) AddMultMV(CONJ(alpha),A.Conjugate(),
	  x.Conjugate(),y.Conjugate());
      else if (A.isrm() && x.step()==1) 
	if (IMAG(alpha) == RealType(T)(0)) 
	  RowMajorAddMultMV(REAL(alpha),A,x,y);
	else RowMajorAddMultMV(alpha,A,x,y);
      else if (A.iscm() && y.step()==1) 
	if (IMAG(alpha) == RealType(T)(0)) 
	  ColMajorAddMultMV(REAL(alpha),A,x,y);
	else ColMajorAddMultMV(alpha,A,x,y);
      else RowAddMultMV(alpha,A,x,y);
    }
#ifdef XDEBUG
    if (Norm(y-y2) > 0.001*Norm(y)) {
      cerr<<"AddMultMV: alpha = "<<alpha<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"x = "<<Type(x)<<" step "<<x.step()<<"  "<<x<<endl;
      cerr<<"y = "<<Type(y)<<" step "<<y.step()<<"  "<<y0<<endl;
      cerr<<"-> y = "<<y<<endl;
      cerr<<"y2 = "<<y2<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta, class Tx> void MultMV(
      const T alpha, const GenUpperTriMatrix<Ta>& A,
      const GenVector<Tx>& x, const T beta, const VectorView<T>& y)
    // y = alpha * A * x + beta * y    
  { 
#ifdef XDEBUG
    Vector<T> y2 = beta*y+alpha*Matrix<Ta>(A)*x;
    Vector<T> y0 = y;
#endif
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());

    if (y.size() > 0) {
      if (alpha==T(0)) {
	y *= beta;
      } else {
	if (beta == T(0)) {
	  y = x;
	  MultEqMV(A,y);
	  y *= alpha;
	} else if (y.SameStorageAs(x)) {
	  Vector<T> temp = x;
	  MultEqMV(A,temp.View());
	  y *= beta;
	  y += alpha * temp;
	} else {
	  y *= beta;
	  AddMultMV(alpha,A,x,y);
	}
      } 
    }
#ifdef XDEBUG
    if (Norm(y-y2) > 0.001*Norm(y)) {
      cerr<<"TriMultMV: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"x = "<<Type(x)<<" step "<<x.step()<<"  "<<x<<endl;
      cerr<<"y = "<<Type(y)<<" step "<<y.step()<<"  "<<y0<<endl;
      cerr<<"-> y = "<<y<<endl;
      cerr<<"y2 = "<<y2<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta, class Tx> void MultMV(
      const T alpha, const GenLowerTriMatrix<Ta>& A,
      const GenVector<Tx>& x, const T beta, const VectorView<T>& y)
    // y = alpha * A * x + beta * y    
  { 
#ifdef XDEBUG
    Vector<T> y2 = beta*y+alpha*Matrix<Ta>(A)*x;
    Vector<T> y0 = y;
#endif

    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());

    if (y.size() > 0) {
      if (alpha==T(0)) {
	y *= beta;
      } else {
	if (beta == T(0)) {
	  y = x;
	  MultEqMV(A,y);
	  y *= alpha;
	} else if (y.SameStorageAs(x)) {
	  Vector<T> temp = x;
	  MultEqMV(A,temp.View());
	  y *= beta;
	  y += alpha * temp;
	} else {
	  y *= beta;
	  AddMultMV(alpha,A,x,y);
	}
      } 
    }
#ifdef XDEBUG
    if (Norm(y-y2) > 0.001*Norm(y)) {
      cerr<<"TriMultMV: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"x = "<<Type(x)<<" step "<<x.step()<<"  "<<x<<endl;
      cerr<<"y = "<<Type(y)<<" step "<<y.step()<<"  "<<y0<<endl;
      cerr<<"-> y = "<<y<<endl;
      cerr<<"y2 = "<<y2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_TriMatrixArith_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


