
#include "TMV_VectorArith_Inline.h"
#include "TMV.h"
#include "TMV_Band.h"

namespace tmv {

  //
  // MultMV
  //

  template <ConjItType Ca, class T1, class Ta, class Tx, ConjItType Cx, class T2, class Ty>
    void RowMajorMultMV(const T1 alpha, 
	const GenBandMatrix<Ta>& A, CVIt<Tx,Unit,Cx> xit,
	const T2 beta, const VectorView<Ty>& y)
    {
      TMVAssert(A.isrm());
      TMVAssert(alpha != T1(0));
      TMVAssert(A.ct() == Ca);

      size_t i=0;
      size_t j1=0;
      size_t j2=A.nhi()+1;
      size_t k=A.nlo();
      size_t len = A.nhi()+1; // len = j2-j1
      const int ds = A.stepi()+1;
      VIt<Ty,Step,NonConj> yit = y.begin();
      const Ta* Aptr = A.cptr();

      if (beta == T2(0)) {
	if (alpha == T1(1)) {
	  for(; i<A.colsize(); ++i, ++yit) {
	    DoMultVV(*yit,CVIt<Ta,Unit,Ca>(Aptr,1),xit,len);
	    if (k>0) {--k; ++len; Aptr+=A.stepi(); }
	    else {++j1; ++xit; Aptr+=ds; }
	    if (j2<A.rowsize()) ++j2;
	    else { --len; if (j1==A.rowsize()) { ++i, ++yit; break; } }
	  }
	}
	else if (alpha == T1(-1)) {
	  for(; i<A.colsize(); ++i, ++yit) {
	    DoMultVV(*yit,CVIt<Ta,Unit,Ca>(Aptr,1),xit,len);
	    *yit = -*yit;
	    if (k>0) {--k; ++len; Aptr+=A.stepi(); }
	    else {++j1; ++xit; Aptr+=ds; }
	    if (j2<A.rowsize()) ++j2;
	    else { --len; if (j1==A.rowsize()) { ++i, ++yit; break; } }
	  }
	}
	else {
	  for(; i<A.colsize(); ++i, ++yit) {
	    DoMultVV(*yit,CVIt<Ta,Unit,Ca>(Aptr,1),xit,len);
	    *yit *= alpha;
	    if (k>0) {--k; ++len; Aptr+=A.stepi(); }
	    else {++j1; ++xit; Aptr+=ds; }
	    if (j2<A.rowsize()) ++j2;
	    else { --len; if (j1==A.rowsize()) { ++i, ++yit; break; } }
	  }
	}
	for(; i<A.colsize(); ++i, ++yit) *yit = Ty(0);
      } else if (beta == T2(1)) {
	if (alpha == T1(1)) {
	  for(; i<A.colsize(); ++i, ++yit) {
	    Ty temp;
	    DoMultVV(temp,CVIt<Ta,Unit,Ca>(Aptr,1),xit,len);
	    *yit += temp;
	    if (k>0) {--k; ++len; Aptr+=A.stepi(); }
	    else {++j1; ++xit; Aptr+=ds; }
	    if (j2<A.rowsize()) ++j2;
	    else { --len; if (j1==A.rowsize()) { ++i, ++yit; break; } }
	  }
	}
	else if (alpha == T1(-1)) {
	  for(; i<A.colsize(); ++i, ++yit) {
	    Ty temp;
	    DoMultVV(temp,CVIt<Ta,Unit,Ca>(Aptr,1),xit,len);
	    *yit -= temp;
	    if (k>0) {--k; ++len; Aptr+=A.stepi(); }
	    else {++j1; ++xit; Aptr+=ds; }
	    if (j2<A.rowsize()) ++j2;
	    else { --len; if (j1==A.rowsize()) { ++i, ++yit; break; } }
	  }
	}
	else {
	  for(; i<A.colsize(); ++i, ++yit) {
	    Ty temp;
	    DoMultVV(temp,CVIt<Ta,Unit,Ca>(Aptr,1),xit,len);
	    *yit += temp*alpha;
	    if (k>0) {--k; ++len; Aptr+=A.stepi(); }
	    else {++j1; ++xit; Aptr+=ds; }
	    if (j2<A.rowsize()) ++j2;
	    else { --len; if (j1==A.rowsize()) { ++i, ++yit; break; } }
	  }
	}
      } else {
	if (alpha == T1(1)) {
	  for(; i<A.colsize(); ++i, ++yit) {
	    Ty temp;
	    DoMultVV(temp,CVIt<Ta,Unit,Ca>(Aptr,1),xit,len);
	    *yit = beta*(*yit) + temp;
	    if (k>0) {--k; ++len; Aptr+=A.stepi(); }
	    else {++j1; ++xit; Aptr+=ds; }
	    if (j2<A.rowsize()) ++j2;
	    else { --len; if (j1==A.rowsize()) { ++i, ++yit; break; } }
	  }
	}
	else if (alpha == T1(-1)) {
	  for(; i<A.colsize(); ++i, ++yit) {
	    Ty temp;
	    DoMultVV(temp,CVIt<Ta,Unit,Ca>(Aptr,1),xit,len);
	    *yit = beta*(*yit) - temp;
	    if (k>0) {--k; ++len; Aptr+=A.stepi(); }
	    else {++j1; ++xit; Aptr+=ds; }
	    if (j2<A.rowsize()) ++j2;
	    else { --len; if (j1==A.rowsize()) { ++i, ++yit; break; } }
	  }
	}
	else {
	  for(; i<A.colsize(); ++i, ++yit) {
	    Ty temp;
	    DoMultVV(temp,CVIt<Ta,Unit,Ca>(Aptr,1),xit,len);
	    *yit = beta*(*yit) + alpha*temp;
	    if (k>0) {--k; ++len; Aptr+=A.stepi(); }
	    else {++j1; ++xit; Aptr+=ds; }
	    if (j2<A.rowsize()) ++j2;
	    else { --len; if (j1==A.rowsize()) { ++i, ++yit; break; } }
	  }
	}
	for(;i<A.colsize(); ++i, ++yit) *yit *= beta;
      }
    }

  template <class T1, class Ta, class Tx, class T2, class Ty> 
    void DoRowMajorMultMV2(const T1 alpha, 
	const GenBandMatrix<Ta>& A, const GenVector<Tx>& x,
	const T2 beta, const VectorView<Ty>& y)
    {
      if (A.isconj())
	if (x.isconj())
	  RowMajorMultMV<Conj>(alpha,A,CVIt<Tx,Unit,Conj>(x.begin()),beta,y);
	else
	  RowMajorMultMV<Conj>(alpha,A,CVIt<Tx,Unit,NonConj>(x.begin()),beta,y);
      else
	if (x.isconj())
	  RowMajorMultMV<NonConj>(alpha,A,CVIt<Tx,Unit,Conj>(x.begin()),beta,y);
	else
	  RowMajorMultMV<NonConj>(alpha,A,CVIt<Tx,Unit,NonConj>(x.begin()),
	      beta,y);
    }

  template <class T, class Ta, class Tx> void DoRowMajorMultMV(const T alpha, 
      const GenBandMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
  {
    if (IMAG(alpha) == RealType(T)(0))
      if (IMAG(beta) == RealType(T)(0))
	DoRowMajorMultMV2(REAL(alpha),A,x,REAL(beta),y);
      else
	DoRowMajorMultMV2(REAL(alpha),A,x,beta,y);
    else
      if (IMAG(beta) == RealType(T)(0))
	DoRowMajorMultMV2(alpha,A,x,REAL(beta),y);
      else
	DoRowMajorMultMV2(alpha,A,x,beta,y);
  }

  template <ConjItType Ca, class T1, class Ta, class Tx, class T2, class Ty>
    void ColMajorMultMV(const T1 alpha, const GenBandMatrix<Ta>& A,
	const GenVector<Tx>& x, const T2 beta, const VectorView<Ty> y)
    {
      TMVAssert(A.iscm());
      TMVAssert(y.step()==1);
      TMVAssert(A.colsize() == y.size());
      TMVAssert(alpha != T1(0));
      TMVAssert(y.size() > 0);
      TMVAssert(y.ct() == NonConj);
      TMVAssert(A.ct() == Ca);

      size_t i1=0;
      size_t i2=A.nlo()+1;
      size_t k=A.nhi();
      size_t len = A.nlo()+1; // = i2-i1
      const int ds = A.stepj()+1;
      VIt<Ty,Unit,NonConj> yit = y.begin();
      CVIter<Tx> xit = x.begin();
      const Ta* Aptr = A.cptr();

      if (beta == T2(0)) {
	y.Zero();
	for(size_t j=0; j<A.rowsize(); ++j,++xit) {
	  if (*xit != Tx(0)) {
	    DoAddVV(*xit,CVIt<Ta,Unit,Ca>(Aptr,1),yit,len);
	  }
	  if (k>0) { --k; Aptr+=A.stepj(); }
	  else { ++i1; ++yit; --len; Aptr+=ds; }
	  if (i2<A.colsize()) {++i2; ++len;}
	  else if (i1==A.colsize()) break;
	}
	if (alpha != T1(1)) y *= alpha;
      } else {
	if (alpha == T1(1)) {
	  if (beta != T2(1)) y *= beta;
	  for(size_t j=0; j<A.rowsize(); ++j,++xit) {
	    if (*xit != Tx(0)) {
	      DoAddVV(*xit,CVIt<Ta,Unit,Ca>(Aptr,1),yit,len);
	    }
	    if (k>0) { --k; Aptr+=A.stepj(); }
	    else { ++i1; ++yit; --len; Aptr+=ds; }
	    if (i2<A.colsize()) {++i2; ++len;}
	    else if (i1==A.colsize()) break;
	  }
	} else if (alpha == T1(-1)) {
	  if (beta != T2(1)) y *= beta;
	  for(size_t j=0; j<A.rowsize(); ++j,++xit) {
	    if (*xit != Tx(0)) {
	      DoAddVV(-*xit,CVIt<Ta,Unit,Ca>(Aptr,1),yit,len);
	    }
	    if (k>0) { --k; Aptr+=A.stepj(); }
	    else { ++i1; ++yit; --len; Aptr+=ds; }
	    if (i2<A.colsize()) {++i2; ++len;}
	    else if (i1==A.colsize()) break;
	  }
	} else {
	  // Requires temporary
	  Vector<Ty> betay = beta*y;
	  ColMajorMultMV<Ca>(alpha,A,x,RealType(Ty)(0),y);
	  y += betay;
	}
      }
    }

  template <class T1, class Ta, class Tx, class T2, class Ty> 
    void DoColMajorMultMV2(const T1 alpha, 
	const GenBandMatrix<Ta>& A, const GenVector<Tx>& x,
	const T2 beta, const VectorView<Ty>& y)
    {
      if (A.isconj()) ColMajorMultMV<Conj>(alpha,A,x,beta,y);
      else ColMajorMultMV<NonConj>(alpha,A,x,beta,y);
    }

  template <class T, class Ta, class Tx> void DoColMajorMultMV(const T alpha, 
      const GenBandMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
  {
    if (IMAG(alpha) == RealType(T)(0))
      if (IMAG(beta) == RealType(T)(0))
	DoColMajorMultMV2(REAL(alpha),A,x,REAL(beta),y);
      else
	DoColMajorMultMV2(REAL(alpha),A,x,beta,y);
    else
      if (IMAG(beta) == RealType(T)(0))
	DoColMajorMultMV2(alpha,A,x,REAL(beta),y);
      else
	DoColMajorMultMV2(alpha,A,x,beta,y);
  }

  template <ConjItType Ca, class T1, class Ta, class Tx, ConjItType Cx, class T2, class Ty>
    void DiagMajorMultMV(const T1 alpha, 
	const GenBandMatrix<Ta>& A, CVIt<Tx,Unit,Cx> xit,
	const T2 beta, const VectorView<Ty>& y)
    {
      TMVAssert(alpha != T1(0));
      TMVAssert(A.colsize() == y.size());
      TMVAssert(A.isdm());
      TMVAssert(y.step()==1);
      TMVAssert(y.size() > 0);
      TMVAssert(A.ct() == Ca);
      TMVAssert(y.ct() == NonConj);

      size_t j2=min(A.colsize()-A.nlo(),A.rowsize());
      size_t len=j2; // == j2-j1
      VIt<Ty,Unit,NonConj> yit = y.begin()+A.nlo();
      const Ta* Aptr = A.cptr()+A.nlo()*A.stepi();
      if (beta == T2(0)) {
	y.Zero();
	for(int k=-A.nlo(); k<=A.nhi(); ++k) {
	  DoAddElementProd(RealType(Ty)(1),CVIt<Ta,Unit,Ca>(Aptr,1),
	      xit,RealType(Ty)(1),yit,len);
	  if (k<0) { --yit; ++len; Aptr-=A.stepi(); } 
	  else { ++xit; Aptr+=A.stepj(); }
	  if (j2 < A.rowsize()) ++j2; else --len; 
	}
	if (alpha != T1(1)) y *= alpha;
      } else {
	if (alpha == T1(1)) {
	  if (beta != T2(1)) y *= beta;
	  for(int k=-A.nlo(); k<=A.nhi(); ++k) {
	    DoAddElementProd(RealType(Ty)(1),CVIt<Ta,Unit,Ca>(Aptr,1),
		xit,RealType(Ty)(1),yit,len);
	    if (k<0) { --yit; ++len; Aptr-=A.stepi(); } 
	    else { ++xit; Aptr+=A.stepj(); }
	    if (j2 < A.rowsize()) ++j2; else --len; 
	  }
	} else if (alpha == T1(-1)) {
	  if (beta != T2(1)) y *= beta;
	  for(int k=-A.nlo(); k<=A.nhi(); ++k) {
	    DoAddElementProd(RealType(Ty)(-1),CVIt<Ta,Unit,Ca>(Aptr,1),
		xit,RealType(Ty)(1),yit,len);
	    if (k<0) { --yit; ++len; Aptr-=A.stepi(); } 
	    else { ++xit; Aptr+=A.stepj(); }
	    if (j2 < A.rowsize()) ++j2; else --len; 
	  }
	} else {
	  // Requires temporary
	  Vector<Ty> betay = beta*y;
	  DiagMajorMultMV<Ca>(alpha,A,xit,RealType(Ty)(0),y);
	  DoAddVV(RealType(Ty)(1),CVIt<Ty,Unit,NonConj>(betay.begin()),
	      VIt<Ty,Unit,NonConj>(y.begin()),y.size());
	}
      }
    }

  template <class T1, class Ta, class Tx, class T2, class Ty> 
    void DoDiagMajorMultMV2(const T1 alpha, 
	const GenBandMatrix<Ta>& A, const GenVector<Tx>& x,
	const T2 beta, const VectorView<Ty>& y)
    {
      if (A.isconj())
	if (x.isconj())
	  DiagMajorMultMV<Conj>(alpha,A,CVIt<Tx,Unit,Conj>(x.begin()),beta,y);
	else
	  DiagMajorMultMV<Conj>(alpha,A,CVIt<Tx,Unit,NonConj>(x.begin()),
	      beta,y);
      else
	if (x.isconj())
	  DiagMajorMultMV<NonConj>(alpha,A,CVIt<Tx,Unit,Conj>(x.begin()),
	      beta,y);
	else
	  DiagMajorMultMV<NonConj>(alpha,A,CVIt<Tx,Unit,NonConj>(x.begin()),
	      beta,y);
    }

  template <class T, class Ta, class Tx> void DoDiagMajorMultMV(const T alpha, 
      const GenBandMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
  {
    if (IMAG(alpha) == RealType(T)(0))
      if (IMAG(beta) == RealType(T)(0))
	DoDiagMajorMultMV2(REAL(alpha),A,x,REAL(beta),y);
      else
	DoDiagMajorMultMV2(REAL(alpha),A,x,beta,y);
    else
      if (IMAG(beta) == RealType(T)(0))
	DoDiagMajorMultMV2(alpha,A,x,REAL(beta),y);
      else
	DoDiagMajorMultMV2(alpha,A,x,beta,y);
  }

  template <class T, class Ta, class Tx> void DiagMultMV(const T alpha, 
      const GenBandMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
  {
    TMVAssert(alpha != T(0));
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(y.ct() == NonConj);

    // A.diag(k) goes from i=i1..i2 and j=j1..j2
    // These start with ranges for k= -A.nlo(r
    size_t j1=0;
    size_t j2=min(A.colsize()-A.nlo(),A.rowsize());
    size_t i1=A.nlo();
    size_t i2=i1+j2;
    if (beta == T(0)) {
      y.Zero();
      for(int k=-A.nlo(); k<=A.nhi(); ++k) {
	AddElementProd(T(1),A.diag(k),x.SubVector(j1,j2),
	    T(1),y.SubVector(i1,i2));
	if (k<0) --i1; else ++j1;
	if (j2 < A.rowsize()) ++j2; else --i2;
      }
      if (alpha != T(1)) y *= alpha;
    } else {
      if (alpha == T(1)) {
	if (beta != T(1)) y *= beta;
	for(int k=-A.nlo(); k<=A.nhi(); ++k) {
	  AddElementProd(T(1),A.diag(k),x.SubVector(j1,j2),
	      T(1),y.SubVector(i1,i2));
	  if (k<0) --i1; else ++j1;
	  if (j2 < A.rowsize()) ++j2; else --i2;
	}
      } else if (alpha == T(-1)) {
	if (beta != T(1)) y *= beta;
	for(int k=-A.nlo(); k<=A.nhi(); ++k) {
	  AddElementProd(T(-1),A.diag(k),x.SubVector(j1,j2),
	      T(1),y.SubVector(i1,i2));
	  if (k<0) --i1; else ++j1;
	  if (j2 < A.rowsize()) ++j2; else --i2;
	}
      } else {
	// Requires temporary
	Vector<T> betay = beta*y;
	DiagMultMV(alpha,A,x,T(0),y);
	y += betay;
      }
    }
  }

  template <class T, class Ta, class Tx> inline void NonBlasMultMV(
      const T alpha, const GenBandMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
    // y = alpha * A * x + beta * y
  {
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"x = "<<Type(x)<<" step "<<x.step()<<"  "<<x<<endl;
    //cerr<<"y = "<<Type(y)<<" step "<<y.step()<<"  "<<y<<endl;
    //Matrix<Ta,RowMajor> A2 = A;
    //Vector<T> y2 = beta*y;
    //y2 += alpha*A2*x;

    TMVAssert(alpha != T(0));
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(y.ct() == NonConj);

    // The main performance issues are to maximize the length
    // of the vectors in the innermost loop and to use
    // unit stride vectors (which are faster to pipeline).
    // It seems that the step size is the more important concern
    // generally, so I check that one first.
    //
    // For BandMatrices, if the strides don't work well, the 
    // Diag algorithm is best, since there are far fewer vector
    // touches in this case.
    //
    // However, as with general Matrices, the column major algorithm,
    // and now also the diag major algorithm requires that either
    // beta == 0 or alpha == +-1.  Otherwise it needs extra storage
    // or extra multiplications. So check for this first.
    //

    if (A.isrm() && x.step()==1) 
      DoRowMajorMultMV(alpha,A,x,beta,y);
    else if (A.iscm() && y.step()==1) 
      DoColMajorMultMV(alpha,A,x,beta,y);
    else if (A.isdm() && x.step()==1 && y.step()==1) 
      DoDiagMajorMultMV(alpha,A,x,beta,y);
    else DiagMultMV(alpha,A,x,beta,y);
    //if (Norm(y2-y) > 0.1) {
      //cerr<<"y = "<<y<<endl;
      //cerr<<"y2 = "<<y2<<endl;
      //abort();
    //}
  }

#ifdef BLAS
  template <class T, class Ta, class Tx> inline void BlasMultMV(
      const T alpha, const GenBandMatrix<Ta>& A,
      const GenVector<Tx>& x, const T beta, const VectorView<T>& y)
  { NonBlasMultMV(alpha,A,x,beta,y); }
  template <> inline void BlasMultMV(const double alpha,
      const GenBandMatrix<double>& A, const GenVector<double>& x,
      const double beta, const VectorView<double>& y)
  {
    TMVAssert(alpha != 0.);
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    if (A.isrm())
      cblas_dgbmv(CblasRowMajor, CblasNoTrans,
	  A.colsize(),A.rowsize(),A.nlo(),A.nhi(),alpha,
	  A.cptr()-A.nlo(), A.stepi()+1,
	  x.cptr(),x.step(),beta,y.ptr(),y.step()); 
    else
      cblas_dgbmv(CblasColMajor, CblasNoTrans,
	  A.colsize(),A.rowsize(),A.nlo(),A.nhi(),alpha,
	  A.cptr()-A.nhi(), A.stepj()+1,
	  x.cptr(),x.step(),beta,y.ptr(),y.step()); 
  }
  template <> inline void BlasMultMV(const complex<double> alpha,
      const GenBandMatrix<complex<double> >& A,
      const GenVector<complex<double> >& x,
      const complex<double> beta, const VectorView<complex<double> >& y)
  {
    TMVAssert(alpha != complex<double>(0));
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(y.ct() == NonConj);
    if (x.isconj()) NonBlasMultMV(alpha,A,x,beta,y);
    else if (A.isconj())
      if (A.isrm()) 
	cblas_zgbmv(CblasColMajor, CblasConjTrans,
	    A.rowsize(),A.colsize(),A.nhi(),A.nlo(),&alpha,
	    A.cptr()-A.nlo(), A.stepi()+1,
	    x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
      else
	cblas_zgbmv(CblasRowMajor, CblasConjTrans,
	    A.rowsize(),A.colsize(),A.nhi(),A.nlo(),&alpha,
	    A.cptr()-A.nhi(), A.stepj()+1,
	    x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
    else
      if (A.isrm())
	cblas_zgbmv(CblasRowMajor, CblasNoTrans,
	    A.colsize(),A.rowsize(),A.nlo(),A.nhi(),&alpha,
	    A.cptr()-A.nlo(), A.stepi()+1,
	    x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
      else
	cblas_zgbmv(CblasColMajor, CblasNoTrans,
	    A.colsize(),A.rowsize(),A.nlo(),A.nhi(),&alpha,
	    A.cptr()-A.nhi(), A.stepj()+1,
	    x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
  }
  template <> inline void BlasMultMV(
      const complex<double> alpha, const GenBandMatrix<double>& A,
      const GenVector<complex<double> >& x,
      const complex<double> beta, const VectorView<complex<double> >& y)
  {
    TMVAssert(y.ct() == NonConj);
    if (IMAG(alpha) == double(0) && IMAG(beta) == double(0)) {
      BlasMultMV(REAL(alpha),A,x.Real(),REAL(beta),y.Real());
      if (x.isconj())
	BlasMultMV(-REAL(alpha),A,x.Conjugate().Imag(),REAL(beta),y.Imag());
      else 
	BlasMultMV(REAL(alpha),A,x.Imag(),REAL(beta),y.Imag());
    } else NonBlasMultMV(alpha,A,x,beta,y);
  }
  template <> inline void BlasMultMV(
      const complex<double> alpha, const GenBandMatrix<double>& A,
      const GenVector<double>& x,
      const complex<double> beta, const VectorView<complex<double> >& y)
  {
    TMVAssert(y.ct() == NonConj);
    if (IMAG(alpha) == double(0) && IMAG(beta) == double(0)) {
      BlasMultMV(REAL(alpha),A,x,REAL(beta),y.Real());
      if (REAL(beta) != double(1)) y.Imag() *= REAL(beta);
    } else NonBlasMultMV(alpha,A,x,beta,y);
  }
#ifndef NOFLOAT
  template <> inline void BlasMultMV(const float alpha,
      const GenBandMatrix<float>& A, const GenVector<float>& x,
      const float beta, const VectorView<float>& y)
  {
    TMVAssert(alpha != 0.);
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    if (A.isrm())
      cblas_sgbmv(CblasRowMajor, CblasNoTrans,
	  A.colsize(),A.rowsize(),A.nlo(),A.nhi(),alpha,
	  A.cptr()-A.nlo(), A.stepi()+1,
	  x.cptr(),x.step(),beta,y.ptr(),y.step()); 
    else
      cblas_sgbmv(CblasColMajor, CblasNoTrans,
	  A.colsize(),A.rowsize(),A.nlo(),A.nhi(),alpha,
	  A.cptr()-A.nhi(), A.stepj()+1,
	  x.cptr(),x.step(),beta,y.ptr(),y.step()); 
  }
  template <> inline void BlasMultMV(const complex<float> alpha,
      const GenBandMatrix<complex<float> >& A,
      const GenVector<complex<float> >& x,
      const complex<float> beta, const VectorView<complex<float> >& y)
  {
    TMVAssert(alpha != complex<float>(0));
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(y.ct() == NonConj);
    if (x.isconj()) NonBlasMultMV(alpha,A,x,beta,y);
    else if (A.isconj())
      if (A.isrm()) 
	cblas_cgbmv(CblasColMajor, CblasConjTrans,
	    A.rowsize(),A.colsize(),A.nhi(),A.nlo(),&alpha,
	    A.cptr()-A.nlo(), A.stepi()+1,
	    x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
      else
	cblas_cgbmv(CblasRowMajor, CblasConjTrans,
	    A.rowsize(),A.colsize(),A.nhi(),A.nlo(),&alpha,
	    A.cptr()-A.nhi(), A.stepj()+1,
	    x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
    else
      if (A.isrm())
	cblas_cgbmv(CblasRowMajor, CblasNoTrans,
	    A.colsize(),A.rowsize(),A.nlo(),A.nhi(),&alpha,
	    A.cptr()-A.nlo(), A.stepi()+1,
	    x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
      else
	cblas_cgbmv(CblasColMajor, CblasNoTrans,
	    A.colsize(),A.rowsize(),A.nlo(),A.nhi(),&alpha,
	    A.cptr()-A.nhi(), A.stepj()+1,
	    x.cptr(),x.step(),&beta,y.ptr(),y.step()); 
  }
  template <> inline void BlasMultMV(
      const complex<float> alpha, const GenBandMatrix<float>& A,
      const GenVector<complex<float> >& x,
      const complex<float> beta, const VectorView<complex<float> >& y)
  {
    TMVAssert(y.ct() == NonConj);
    if (IMAG(alpha) == float(0) && IMAG(beta) == float(0)) {
      BlasMultMV(REAL(alpha),A,x.Real(),REAL(beta),y.Real());
      if (x.isconj())
	BlasMultMV(-REAL(alpha),A,x.Conjugate().Imag(),REAL(beta),y.Imag());
      else 
	BlasMultMV(REAL(alpha),A,x.Imag(),REAL(beta),y.Imag());
    } else NonBlasMultMV(alpha,A,x,beta,y);
  }
  template <> inline void BlasMultMV(
      const complex<float> alpha, const GenBandMatrix<float>& A,
      const GenVector<float>& x,
      const complex<float> beta, const VectorView<complex<float> >& y)
  {
    TMVAssert(y.ct() == NonConj);
    if (IMAG(alpha) == float(0) && IMAG(beta) == float(0)) {
      BlasMultMV(REAL(alpha),A,x,REAL(beta),y.Real());
      if (REAL(beta) != float(1)) y.Imag() *= REAL(beta);
    } else NonBlasMultMV(alpha,A,x,beta,y);
  }
#endif // NOFLOAT
#endif // BLAS

  template <class T, class Ta, class Tx> inline void DoMultMV(const T alpha,
      const GenBandMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(x.size() > 0);
    TMVAssert(y.size() > 0);

    if (y.isconj()) DoMultMV(CONJ(alpha),A.QuickConjugate(),
	x.Conjugate(),CONJ(beta),y.Conjugate());
    else {
#ifdef BLAS
      if ((A.isrm() || A.iscm()) && x.step()>0 && y.step()>0) {
	if (A.rowsize() > A.colsize()) {
	  size_t N1 = A.colsize() + A.nhi();
	  if (A.rowsize() > N1) {
	    BlasMultMV(alpha,A.SubBandMatrix(0,A.colsize(),0,N1,
		  A.nlo(),A.nhi()),x.SubVector(0,N1),beta,y);
	  }
	  else BlasMultMV(alpha,A,x,beta,y);
	}
	else if (A.colsize() > A.rowsize()) {
	  size_t N2 = A.rowsize() + A.nlo();
	  if (A.colsize() > N2) {
	    BlasMultMV(alpha,A.SubBandMatrix(0,N2,0,A.rowsize(),
		  A.nlo(),A.nhi()),x,beta,y.SubVector(0,N2));
	    y.SubVector(N2,y.size()).Zero();
	  }
	  else BlasMultMV(alpha,A,x,beta,y);
	}
	else BlasMultMV(alpha,A,x,beta,y);
      }
      else
#endif
	NonBlasMultMV(alpha,A,x,beta,y);
    }
  }

  template <class T, class Ta, class Tx> inline void MultMV(const T alpha,
      const GenBandMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
    // y = alpha * A * x + beta * y    
  { 
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());
    //cerr<<"Start Band: MultMV\n";
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"x = "<<Type(x)<<"  step "<<x.step()<<"  "<<x<<endl;
    //cerr<<"y = "<<Type(y)<<"  step "<<y.step()<<"  "<<y<<endl;
    //Vector<T> y2 = y;
    //y2 *= beta;
    //y2 += alpha*Matrix<T>(A)*x;

    if (y.size() > 0) {
      if (x.size()==0 || alpha==T(0)) MultXV(beta,y);
      else if (x.SameStorageAs(y)) {
	Vector<Tx> temp = x;
	DoMultMV(alpha,A,temp,beta,y);
      } else {
	DoMultMV(alpha,A,x,beta,y);
      } 
    }
    //cerr<<"->y = "<<y<<endl;
    //cerr<<"y2 = "<<y2<<endl;
  }

#define InstFile "TMV_BandMatrixArith_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


