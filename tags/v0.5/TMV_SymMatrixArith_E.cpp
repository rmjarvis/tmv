
#include "TMV_VectorArith_Inline.h"
#include "TMV.h"
#include "TMV_Sym.h"

//#define XDEBUG

namespace tmv {

  // 
  // Rank1Update
  //

  template <class T, class T1, class Tx> void RowMajorRank1Update(
      const T1 alpha, const GenVector<Tx>& x, const SymMatrixView<T>& A)
  {
#ifdef XDEBUG
    //cerr<<"RowMajor Rank1\n";
#endif
    TMVAssert(A.isrm());
    TMVAssert(x.step()==1);
    TMVAssert(A.size() == x.size());
    TMVAssert(alpha != T1(0));
    TMVAssert(IsReal(T1()) || !A.isherm());
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);

    CVIter<Tx> xit = x.begin();
    T* Aptr = A.ptr();
    for (size_t i=0;i<A.size();++i,++xit,Aptr+=A.stepi()) {
      if (*xit != Tx(0)) {
	T ax = *xit;
	if (alpha != T1(1)) ax *= alpha;
	// A.row(i,0,i+1) += ax * A.isherm() ?
	//      x.SubVector(0,i+1).Conjugate() : x.SubVector(0,i+1);;
	if (IMAG(ax) == RealType(T)(0))
	  if (x.isconj() != A.isherm())
	    DoAddVV(REAL(ax),CVIt<Tx,Unit,Conj>(x.cptr(),1),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),i+1);
	  else
	    DoAddVV(REAL(ax),CVIt<Tx,Unit,NonConj>(x.cptr(),1),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),i+1);
	else
	  if (x.isconj() != A.isherm())
	    DoAddVV(ax,CVIt<Tx,Unit,Conj>(x.cptr(),1),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),i+1);
	  else
	    DoAddVV(ax,CVIt<Tx,Unit,NonConj>(x.cptr(),1),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),i+1);
      }
    }
  }

  template <class T, class Tx> inline void DoRowMajorRank1Update(
      const T alpha, const GenVector<Tx>& x, const SymMatrixView<T>& A)
  {
    if (IMAG(alpha) == RealType(T)(0))
      RowMajorRank1Update(REAL(alpha),x,A);
    else
      RowMajorRank1Update(alpha,x,A);
  }

  template <class T, class Tx> void RowRank1Update(
      const T alpha, const GenVector<Tx>& x, const SymMatrixView<T>& A)
  {
#ifdef XDEBUG
    //cerr<<"Row Rank1\n";
#endif
    TMVAssert(A.size() == x.size());
    TMVAssert(alpha != T(0));
    TMVAssert(IMAG(alpha)==RealType(T)(0) || !A.isherm());
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);

    CVIter<Tx> xit = x.begin();
    for (size_t i=0;i<A.size();++i,++xit) {
      T ax = *xit;
      if (alpha != T(1)) ax *= alpha;
      if (A.isherm())
	A.row(i,0,i+1) += ax * x.SubVector(0,i+1).Conjugate();
      else
	A.row(i,0,i+1) += ax * x.SubVector(0,i+1);
    }
  }

  template <class T, class T1, class Tx> void ColMajorRank1Update(
      const T1 alpha, const GenVector<Tx>& x, const SymMatrixView<T>& A)
  {
#ifdef XDEBUG
    //cerr<<"ColMajor Rank1\n";
#endif
    TMVAssert(A.iscm());
    TMVAssert(x.step()==1);
    TMVAssert(A.size() == x.size());
    TMVAssert(alpha != T1(0));
    TMVAssert(IMAG(alpha)==RealType(T1)(0) || !A.isherm());
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);

    CVIter<Tx> xit = x.begin();
    T* Aptr = A.ptr();
    size_t len = A.size();
    const size_t ds = A.stepj()+1;
    for (size_t j=0;j<A.size();++j,++xit,--len,Aptr+=ds) {
      if (*xit != Tx(0)) {
	T ax = A.isherm() ? CONJ(*xit) : *xit;
	if (alpha != T1(1)) ax *= alpha;
	// A.col(j,j,N) += ax * x.SubVector(j,N);
	if (IMAG(ax) == RealType(T)(0))
	  if (x.isconj())
	    DoAddVV(REAL(ax),CVIt<Tx,Unit,Conj>(xit),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),len);
	  else
	    DoAddVV(REAL(ax),CVIt<Tx,Unit,NonConj>(xit),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),len);
	else
	  if (x.isconj())
	    DoAddVV(ax,CVIt<Tx,Unit,Conj>(xit),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),len);
	  else
	    DoAddVV(ax,CVIt<Tx,Unit,NonConj>(xit),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),len);
      }
    }
  }

  template <class T, class Tx> inline void DoColMajorRank1Update(
      const T alpha, const GenVector<Tx>& x, const SymMatrixView<T>& A)
  {
    if (IMAG(alpha) == RealType(T)(0))
      ColMajorRank1Update(REAL(alpha),x,A);
    else
      ColMajorRank1Update(alpha,x,A);
  }

  template <class T, class Tx> void ColRank1Update(
      const T alpha, const GenVector<Tx>& x, const SymMatrixView<T>& A)
  {
#ifdef XDEBUG
    //cerr<<"Col Rank1\n";
#endif
    TMVAssert(A.size() == x.size());
    TMVAssert(alpha != T(0));
    TMVAssert(IMAG(alpha)==RealType(T)(0) || !A.isherm());
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);

    const size_t N = A.size();

    CVIter<Tx> xit = x.begin();
    for (size_t j=0;j<N;++j,++xit) {
      T ax = A.isherm() ? CONJ(*xit) : *xit;
      if (alpha != T(1)) ax *= alpha;
      A.col(j,j,N) += ax * x.SubVector(j,N);
    }
  }

  template <class T, class Tx> void NonBlasRank1Update(
      const T alpha, const GenVector<Tx>& x, const SymMatrixView<T>& A)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(alpha != T(0));
    TMVAssert(IMAG(alpha)==RealType(T)(0) || !A.isherm());
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);

    if (x.step()==1)
      if (A.isrm()) DoRowMajorRank1Update(alpha,x,A);
      else if (A.iscm()) DoColMajorRank1Update(alpha,x,A);
      else ColRank1Update(alpha,x,A);
    else if (A.isrm()) RowRank1Update(alpha,x,A);
    else ColRank1Update(alpha,x,A);
  }

#ifdef BLAS
  template <class T, class Tx> inline void BlasRank1Update(
      const T alpha, const GenVector<Tx>& x, const SymMatrixView<T>& A)
  { NonBlasRank1Update(alpha,x,A); }
  template <> inline void BlasRank1Update(
      const double alpha, const GenVector<double>& x,
      const SymMatrixView<double>& A)
  {
#ifdef XDEBUG
    //cerr<<"Blas Rank1 double\n";
#endif
    TMVAssert(A.size() == x.size());
    TMVAssert(alpha != double(0));
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.sym() == Sym);
    TMVAssert(x.step() > 0);

    if (A.isrm())
      cblas_dsyr(CblasRowMajor,CblasLower,A.size(),alpha,
	  x.cptr(),x.step(),A.ptr(),A.stepi()); 
    else
      cblas_dsyr(CblasColMajor,CblasLower,A.size(),alpha,
	  x.cptr(),x.step(),A.ptr(),A.stepj()); 
  }
  template <> void BlasRank1Update(
      const complex<double> alpha, const GenVector<complex<double> >& x, 
      const SymMatrixView<complex<double> >& A)
  {
#ifdef XDEBUG
    //cerr<<"Blas Rank1 c double\n";
#endif
    TMVAssert(A.size() == x.size());
    TMVAssert(alpha != double(0));
    TMVAssert(IMAG(alpha)==double(0) || !A.isherm());
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(x.step() > 0);

    if (A.isherm()) {
      TMVAssert(IMAG(alpha)==double(0));
      if (A.isrm())
	if (x.isconj()) 
	  if (x.step() == 1) {
	    cblas_zherk(CblasRowMajor,CblasLower,CblasConjTrans,
		A.size(),1,REAL(alpha),x.cptr(),x.size(),
		double(1),A.ptr(),A.stepi()); 
	  }
	  else NonBlasRank1Update(alpha,x,A);
	else 
	  cblas_zher(CblasRowMajor,CblasLower,A.size(),REAL(alpha),
	      x.cptr(),x.step(),A.ptr(),A.stepi()); 
      else
	if (x.isconj()) NonBlasRank1Update(alpha,x,A);
	else 
	  cblas_zher(CblasColMajor,CblasLower,A.size(),REAL(alpha),
	      x.cptr(),x.step(),A.ptr(),A.stepj()); 
    } else {
      if (x.isconj()) NonBlasRank1Update(alpha,x,A);
      else if (x.step() == 1) {
	complex<double> beta(1);
	if (A.isrm())
	  cblas_zsyrk(CblasRowMajor,CblasLower,CblasTrans,
	      A.size(),1,&alpha,x.cptr(),x.size(),&beta,A.ptr(),A.stepi()); 
	else
	  cblas_zsyrk(CblasColMajor,CblasLower,CblasNoTrans,
	      A.size(),1,&alpha,x.cptr(),x.size(),&beta,A.ptr(),A.stepj()); 
      }
      else NonBlasRank1Update(alpha,x,A);
    }
  }
#ifndef NOFLOAT
  template <> inline void BlasRank1Update(
      const float alpha, const GenVector<float>& x,
      const SymMatrixView<float>& A)
  {
#ifdef XDEBUG
    //cerr<<"Blas Rank1 float\n";
#endif
    TMVAssert(A.size() == x.size());
    TMVAssert(alpha != float(0));
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.sym() == Sym);
    TMVAssert(x.step() > 0);

    if (A.isrm())
      cblas_ssyr(CblasRowMajor,CblasLower,A.size(),alpha,
	  x.cptr(),x.step(),A.ptr(),A.stepi()); 
    else
      cblas_ssyr(CblasColMajor,CblasLower,A.size(),alpha,
	  x.cptr(),x.step(),A.ptr(),A.stepj()); 
  }
  template <> void BlasRank1Update(
      const complex<float> alpha, const GenVector<complex<float> >& x, 
      const SymMatrixView<complex<float> >& A)
  {
#ifdef XDEBUG
    //cerr<<"Blas Rank1 c float\n";
#endif
    TMVAssert(A.size() == x.size());
    TMVAssert(alpha != float(0));
    TMVAssert(IMAG(alpha)==float(0) || !A.isherm());
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(x.step() > 0);

    if (A.isherm()) {
      TMVAssert(IMAG(alpha)==float(0));
      if (A.isrm())
	if (x.isconj()) 
	  if (x.step() == 1) {
	    cblas_cherk(CblasRowMajor,CblasLower,CblasConjTrans,
		A.size(),1,REAL(alpha),x.cptr(),x.size(),
		float(1),A.ptr(),A.stepi()); 
	  }
	  else NonBlasRank1Update(alpha,x,A);
	else 
	  cblas_cher(CblasRowMajor,CblasLower,A.size(),REAL(alpha),
	      x.cptr(),x.step(),A.ptr(),A.stepi()); 
      else
	if (x.isconj()) NonBlasRank1Update(alpha,x,A);
	else 
	  cblas_cher(CblasColMajor,CblasLower,A.size(),REAL(alpha),
	      x.cptr(),x.step(),A.ptr(),A.stepj()); 
    } else {
      if (x.isconj()) NonBlasRank1Update(alpha,x,A);
      else if (x.step() == 1) {
	complex<float> beta(1);
	if (A.isrm())
	  cblas_csyrk(CblasRowMajor,CblasLower,CblasTrans,
	      A.size(),1,&alpha,x.cptr(),x.size(),&beta,A.ptr(),A.stepi()); 
	else
	  cblas_csyrk(CblasColMajor,CblasLower,CblasNoTrans,
	      A.size(),1,&alpha,x.cptr(),x.size(),&beta,A.ptr(),A.stepj()); 
      }
      else NonBlasRank1Update(alpha,x,A);
    }
  }
#endif // NOFLOAT
#endif // BLAS

  template <class T, class Tx> void Rank1Update(
      const T alpha, const GenVector<Tx>& x, const SymMatrixView<T>& A)
    // A = A + alpha * x * xT
  {
#ifdef XDEBUG
    Matrix<T> A2 = Matrix<T>(A);
    if (A.isherm())
      A2 += (alpha*x^x.Conjugate());
    else 
      A2 += (alpha*x^x);
    Matrix<T> A0 = A;
#endif

    TMVAssert(A.size() == x.size());
    TMVAssert(IMAG(alpha)==RealType(T)(0) || !A.isherm());
    if (alpha != T(0) && A.size() > 0) {
      if (A.isconj()) 
	Rank1Update(CONJ(alpha),x.Conjugate(),A.Conjugate());
      else if (A.uplo() == Upper) {
	if (A.isherm()) Rank1Update(alpha,x,A.Adjoint());
	else Rank1Update(alpha,x,A.Transpose());
      }
#ifdef BLAS
      else if ((A.isrm() || A.iscm()) && x.step()>0)
	BlasRank1Update(alpha,x,A);
#endif
      else NonBlasRank1Update(alpha,x,A);
    }
  
#ifdef XDEBUG
    if (Norm(A2-A) > 0.001*Norm(A)) {
      cerr<<"Rank1Update: alpha = "<<alpha<<endl;
      cerr<<"x = "<<Type(x)<<"  step = "<<x.step()<<"  "<<x<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> A = "<<A<<endl;
      cerr<<"A2 = "<<A2<<endl;
      abort();
    }
#endif
  }

  // 
  // Rank2Update
  //

  template <class T, class T1, class Tx, class Ty> void RowMajorRank2Update(
      const T1 alpha, const GenVector<Tx>& x, const GenVector<Ty>& y,
      const SymMatrixView<T>& A)
  {
#ifdef XDEBUG
    //cerr<<"RowMajor Rank2\n";
#endif
    TMVAssert(A.isrm());
    TMVAssert(x.step()==1);
    TMVAssert(y.step()==1);
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != T1(0));
    TMVAssert(A.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);

    CVIter<Tx> xit = x.begin();
    CVIter<Ty> yit = y.begin();
    T* Aptr = A.ptr();
    for (size_t i=0;i<A.size();++i,++xit,++yit,Aptr+=A.stepi()) {
      if (*xit != Tx(0)) {
	T ax = *xit;
	if (alpha != T1(1)) ax *= alpha;
	// A.row(i,0,i+1) += ax * A.isherm() ?
	//      y.SubVector(0,i+1).Conjugate() : y.SubVector(0,i+1);
	if (IMAG(ax) == RealType(T)(0))
	  if (IsComplex(Ty()) && y.isconj() != A.isherm())
	    DoAddVV(REAL(ax),CVIt<Ty,Unit,Conj>(y.begin()),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),i+1);
	  else
	    DoAddVV(REAL(ax),CVIt<Ty,Unit,NonConj>(y.begin()),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),i+1);
	else
	  if (IsComplex(Ty()) && y.isconj() != A.isherm())
	    DoAddVV(ax,CVIt<Ty,Unit,Conj>(y.begin()),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),i+1);
	  else
	    DoAddVV(ax,CVIt<Ty,Unit,NonConj>(y.begin()),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),i+1);
      }
      if (*yit != Ty(0)) {
	T ay = *yit;
	if (alpha != T1(1)) ay *= A.isherm() ? CONJ(alpha) : alpha;
	// A.row(i,0,i+1) += ay * A.isherm() ?
	//      x.SubVector(0,i+1).Conjugate() : x.SubVector(0,i+1);
	if (IMAG(ay) == RealType(T)(0))
	  if (IsComplex(Tx()) && x.isconj() != A.isherm())
	    DoAddVV(REAL(ay),CVIt<Tx,Unit,Conj>(x.begin()),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),i+1);
	  else
	    DoAddVV(REAL(ay),CVIt<Tx,Unit,NonConj>(x.begin()),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),i+1);
	else
	  if (IsComplex(Tx()) && x.isconj() != A.isherm())
	    DoAddVV(ay,CVIt<Tx,Unit,Conj>(x.begin()),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),i+1);
	  else
	    DoAddVV(ay,CVIt<Tx,Unit,NonConj>(x.begin()),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),i+1);
      }
    }
  }

  template <class T, class Tx, class Ty> inline void DoRowMajorRank2Update(
      const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y,
      const SymMatrixView<T>& A)
  {
    if (IMAG(alpha) == RealType(T)(0))
      RowMajorRank2Update(REAL(alpha),x,y,A);
    else
      RowMajorRank2Update(alpha,x,y,A);
  }

  template <class T, class Tx, class Ty> void RowRank2Update(
      const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y,
      const SymMatrixView<T>& A)
  {
#ifdef XDEBUG
    //cerr<<"Row Rank2\n";
#endif
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);

    CVIter<Tx> xit = x.begin();
    CVIter<Ty> yit = y.begin();
    for (size_t i=0;i<A.size();++i,++xit,++yit) {
      T ax = *xit;
      if (alpha != T(1)) ax *= alpha;
      if (A.isherm())
	A.row(i,0,i+1) += ax * y.SubVector(0,i+1).Conjugate();
      else
	A.row(i,0,i+1) += ax * y.SubVector(0,i+1);
      T ay = *yit;
      if (alpha != T(1)) ay *= A.isherm() ? CONJ(alpha) : alpha;
      if (A.isherm())
	A.row(i,0,i+1) += ay * x.SubVector(0,i+1).Conjugate();
      else
	A.row(i,0,i+1) += ay * x.SubVector(0,i+1);
    }
  }

  template <class T, class T1, class Tx, class Ty> void ColMajorRank2Update(
      const T1 alpha, const GenVector<Tx>& x, const GenVector<Ty>& y,
      const SymMatrixView<T>& A)
  {
#ifdef XDEBUG
    //cerr<<"ColMajor Rank2\n";
#endif
    TMVAssert(A.iscm());
    TMVAssert(x.step()==1);
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != T1(0));
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);

    CVIter<Tx> xit = x.begin();
    CVIter<Ty> yit = y.begin();
    T* Aptr = A.ptr();
    size_t len = A.size();
    const size_t ds = A.stepj()+1;
    for (size_t j=0;j<A.size();++j,++xit,++yit,--len,Aptr+=ds) {
      if (*yit != Ty(0)) {
	T ay = A.isherm() ? CONJ(*yit) : *yit;
	if (alpha != T1(1)) ay *= alpha;
	// A.col(j,j,N) += ay * x.SubVector(j,N);
	if (IMAG(ay) == RealType(T)(0))
	  if (x.isconj())
	    DoAddVV(REAL(ay),CVIt<Tx,Unit,Conj>(xit),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),len);
	  else
	    DoAddVV(REAL(ay),CVIt<Tx,Unit,NonConj>(xit),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),len);
	else
	  if (x.isconj())
	    DoAddVV(ay,CVIt<Tx,Unit,Conj>(xit),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),len);
	  else
	    DoAddVV(ay,CVIt<Tx,Unit,NonConj>(xit),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),len);
      }
      if (*xit != Tx(0)) {
	T ax = A.isherm() ? CONJ(*xit) : *xit;
	if (alpha != T1(1)) ax *= A.isherm() ? CONJ(alpha) : alpha;
	// A.col(j,j,N) += ax * y.SubVector(j,N);
	if (IMAG(ax) == RealType(T)(0))
	  if (x.isconj())
	    DoAddVV(REAL(ax),CVIt<Ty,Unit,Conj>(yit),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),len);
	  else
	    DoAddVV(REAL(ax),CVIt<Ty,Unit,NonConj>(yit),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),len);
	else
	  if (x.isconj())
	    DoAddVV(ax,CVIt<Ty,Unit,Conj>(yit),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),len);
	  else
	    DoAddVV(ax,CVIt<Ty,Unit,NonConj>(yit),
		VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),len);
      }
    }
  }

  template <class T, class Tx, class Ty> inline void DoColMajorRank2Update(
      const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y,
      const SymMatrixView<T>& A)
  {
    if (IMAG(alpha) == RealType(T)(0))
      ColMajorRank2Update(REAL(alpha),x,y,A);
    else
      ColMajorRank2Update(alpha,x,y,A);
  }

  template <class T, class Tx, class Ty> void ColRank2Update(
      const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y,
      const SymMatrixView<T>& A)
  {
#ifdef XDEBUG
    //cerr<<"Col Rank2\n";
#endif
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);

    const size_t N = A.size();

    CVIter<Tx> xit = x.begin();
    CVIter<Ty> yit = y.begin();
    for (size_t j=0;j<N;++j,++xit) {
      T ay = A.isherm() ? CONJ(*yit) : *yit;
      if (alpha != T(1)) ay *= alpha;
      A.col(j,j,N) += ay * x.SubVector(j,N);
      T ax = A.isherm() ? CONJ(*xit) : *xit;
      if (alpha != T(1)) ax *= A.isherm() ? CONJ(alpha) : alpha;
      A.col(j,j,N) += ax * y.SubVector(j,N);
    }
  }

  template <class T, class Tx, class Ty> void NonBlasRank2Update(
      const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y,
      const SymMatrixView<T>& A)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != T(0));
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);

    if (x.step()==1 && y.step()==1)
      if (A.isrm()) DoRowMajorRank2Update(alpha,x,y,A);
      else if (A.iscm()) DoColMajorRank2Update(alpha,x,y,A);
      else ColRank2Update(alpha,x,y,A);
    else if (A.isrm()) RowRank2Update(alpha,x,y,A);
    else ColRank2Update(alpha,x,y,A);
  }

#ifdef BLAS
  template <class T, class Tx, class Ty> inline void BlasRank2Update(
      const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y,
      const SymMatrixView<T>& A)
  { NonBlasRank2Update(alpha,x,y,A); }
  template <> inline void BlasRank2Update(
      const double alpha, const GenVector<double>& x,
      const GenVector<double>& y, const SymMatrixView<double>& A)
  {
#ifdef XDEBUG
    //cerr<<"Blas Rank2 double\n";
#endif
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != double(0));
    TMVAssert(x.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.sym() == Sym);

    if (A.isrm())
      cblas_dsyr2(CblasRowMajor,CblasLower,A.size(),alpha,
	  x.cptr(),x.step(),y.cptr(),y.step(),A.ptr(),A.stepi()); 
    else
      cblas_dsyr2(CblasColMajor,CblasLower,A.size(),alpha,
	  x.cptr(),x.step(),y.cptr(),y.step(),A.ptr(),A.stepj()); 
  }
  template <> void BlasRank2Update(
      const complex<double> alpha, const GenVector<complex<double> >& x, 
      const GenVector<double>& y, const SymMatrixView<complex<double> >& A)
  {
#ifdef XDEBUG
    //cerr<<"Blas Rank2 c double\n";
#endif
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != double(0));
    TMVAssert(x.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);

    if (x.isconj() != y.isconj())
      NonBlasRank2Update(alpha,x,y,A);
    else if (A.isherm()) {
      if (A.isrm())
	if (x.isconj()) 
	  if (x.step() == 1 && y.step() == 1) {
	    cblas_zher2k(CblasRowMajor,CblasLower,CblasConjTrans,
		A.size(),1,&alpha,y.cptr(),y.size(),x.cptr(),x.size(),
		double(1),A.ptr(),A.stepi()); 
	  }
	  else NonBlasRank2Update(alpha,x,y,A);
	else 
	  cblas_zher2(CblasRowMajor,CblasLower,A.size(),&alpha,
	      x.cptr(),x.step(),y.cptr(),y.step(),A.ptr(),A.stepi()); 
      else
	if (x.isconj()) NonBlasRank2Update(alpha,x,y,A);
	else 
	  cblas_zher2(CblasColMajor,CblasLower,A.size(),&alpha,
	      x.cptr(),x.step(),y.cptr(),y.step(),A.ptr(),A.stepj()); 
    } else {
      if (x.isconj()) NonBlasRank2Update(alpha,x,y,A);
      else if (x.step() == 1 && y.step() == 1) {
	complex<double> beta(1);
	if (A.isrm())
	  cblas_zsyr2k(CblasRowMajor,CblasLower,CblasTrans,
	      A.size(),1,&alpha,x.cptr(),x.size(),y.cptr(),y.size(),
	      &beta,A.ptr(),A.stepi()); 
	else
	  cblas_zsyr2k(CblasColMajor,CblasLower,CblasNoTrans,
	      A.size(),1,&alpha,x.cptr(),x.size(),y.cptr(),y.size(),
	      &beta,A.ptr(),A.stepj()); 
      }
      else NonBlasRank2Update(alpha,x,y,A);
    }
  }
#ifndef NOFLOAT
  template <> inline void BlasRank2Update(
      const float alpha, const GenVector<float>& x,
      const GenVector<float>& y, const SymMatrixView<float>& A)
  {
#ifdef XDEBUG
    //cerr<<"Blas Rank2 float\n";
#endif
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != float(0));
    TMVAssert(x.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);

    if (A.isrm())
      cblas_ssyr2(CblasRowMajor,CblasLower,A.size(),alpha,
	  x.cptr(),x.step(),y.cptr(),y.step(),A.ptr(),A.stepi()); 
    else
      cblas_ssyr2(CblasColMajor,CblasLower,A.size(),alpha,
	  x.cptr(),x.step(),y.cptr(),y.step(),A.ptr(),A.stepj()); 
  }
  template <> void BlasRank2Update(
      const complex<float> alpha, const GenVector<complex<float> >& x, 
      const GenVector<float>& y, const SymMatrixView<complex<float> >& A)
  {
#ifdef XDEBUG
    //cerr<<"Blas Rank2 c float\n";
#endif
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(alpha != float(0));
    TMVAssert(x.size() > 0);
    TMVAssert(x.step() > 0);
    TMVAssert(y.step() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);

    if (x.isconj() != y.isconj())
      NonBlasRank2Update(alpha,x,y,A);
    else if (A.isherm()) {
      if (A.isrm())
	if (x.isconj()) 
	  if (x.step() == 1 && y.step() == 1) {
	    cblas_cher2k(CblasRowMajor,CblasLower,CblasConjTrans,
		A.size(),1,&alpha,y.cptr(),y.size(),x.cptr(),x.size(),
		float(1),A.ptr(),A.stepi()); 
	  }
	  else NonBlasRank2Update(alpha,x,y,A);
	else 
	  cblas_cher2(CblasRowMajor,CblasLower,A.size(),&alpha,
	      x.cptr(),x.step(),y.cptr(),y.step(),A.ptr(),A.stepi()); 
      else
	if (x.isconj()) NonBlasRank2Update(alpha,x,y,A);
	else 
	  cblas_cher2(CblasColMajor,CblasLower,A.size(),&alpha,
	      x.cptr(),x.step(),y.cptr(),y.step(),A.ptr(),A.stepj()); 
    } else {
      if (x.isconj()) NonBlasRank2Update(alpha,x,y,A);
      else if (x.step() == 1 && y.step() == 1) {
	complex<float> beta(1);
	if (A.isrm())
	  cblas_csyr2k(CblasRowMajor,CblasLower,CblasTrans,
	      A.size(),1,&alpha,x.cptr(),x.size(),y.cptr(),y.size(),
	      &beta,A.ptr(),A.stepi()); 
	else
	  cblas_csyr2k(CblasColMajor,CblasLower,CblasNoTrans,
	      A.size(),1,&alpha,x.cptr(),x.size(),y.cptr(),y.size(),
	      &beta,A.ptr(),A.stepj()); 
      }
      else NonBlasRank2Update(alpha,x,y,A);
    }
  }
#endif // NOFLOAT
#endif // BLAS

  template <class T, class Tx, class Ty> void Rank2Update(
      const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y,
      const SymMatrixView<T>& A)
    // if A is sym:  A = A + alpha * (x ^ y + y ^ x)
    // if A is herm: A = A + alpha * x ^ y* + conj(alpha) * y ^ x*
  {
#ifdef XDEBUG
    Matrix<T> A2 = Matrix<T>(A);
    if (A.isherm())
      A2 += (alpha*x^y.Conjugate())+(CONJ(alpha)*y^x.Conjugate());
    else
      A2 += alpha*((x^y) + (y^x));
    Matrix<T> A0 = A;
#endif

    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    if (alpha != T(0) && A.size() > 0) {
      if (A.isconj()) 
	Rank2Update(CONJ(alpha),x.Conjugate(),y.Conjugate(),
	    A.Conjugate());
      else if (A.uplo() == Upper) {
	if (A.isherm()) Rank2Update(alpha,x,y,A.Adjoint());
	else Rank2Update(alpha,x,y,A.Transpose());
      }
#ifdef BLAS
      else if ((A.isrm() || A.iscm()) && x.step()>0 && y.step()>0)
	BlasRank2Update(alpha,x,y,A);
#endif
      else NonBlasRank2Update(alpha,x,y,A);
    }
  
#ifdef XDEBUG
    if (Norm(A2-A) > 0.001*Norm(A)) {
      cerr<<"Rank2Update: alpha = "<<alpha<<endl;
      cerr<<"x = "<<Type(x)<<"  step = "<<x.step()<<"  "<<x<<endl;
      cerr<<"y = "<<Type(y)<<"  step = "<<y.step()<<"  "<<y<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> A = "<<A<<endl;
      cerr<<"A2 = "<<A2<<endl;
      abort();
    }
#endif
  }

  // 
  // RankKUpdate
  //

  template <class T, class Tx> void RowRankKUpdate(
      const T alpha, const GenMatrix<Tx>& x, const SymMatrixView<T>& A)
  {
#ifdef XDEBUG
    cerr<<"Row RankK\n";
#endif
    TMVAssert(x.isrm());
    TMVAssert(A.size() == x.colsize());
    TMVAssert(alpha != T(0));
    TMVAssert(IMAG(alpha)==RealType(T)(0) || !A.isherm());
    TMVAssert(x.colsize() > 0);
    TMVAssert(x.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);

    for (size_t i=0;i<A.size();++i) {
      if (A.isherm())
	A.row(i,0,i+1) += alpha * x.row(i) * x.Rows(0,i+1).Adjoint();
      else
	A.row(i,0,i+1) += alpha * x.row(i) * x.Rows(0,i+1).Transpose();
    }
  }

  template <class T, class Tx> void ColRankKUpdate(
      const T alpha, const GenMatrix<Tx>& x, const SymMatrixView<T>& A)
  {
#ifdef XDEBUG
    cerr<<"Col RankK\n";
#endif
    TMVAssert(x.isrm());
    TMVAssert(A.size() == x.colsize());
    TMVAssert(alpha != T(0));
    TMVAssert(IMAG(alpha)==RealType(T)(0) || !A.isherm());
    TMVAssert(x.colsize() > 0);
    TMVAssert(x.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);

    const size_t N = A.size();

    for (size_t j=0;j<N;++j) {
      if (A.isherm())
	A.col(j,j,N) += alpha * x.Rows(j,N) * x.row(j).Conjugate();
      else
	A.col(j,j,N) += alpha * x.Rows(j,N) * x.row(j);
    }
  }

  template <class T, class Tx> void NonBlasRankKUpdate(
      const T alpha, const GenMatrix<Tx>& x, const SymMatrixView<T>& A)
  {
    TMVAssert(A.size() == x.colsize());
    TMVAssert(alpha != T(0));
    TMVAssert(IMAG(alpha)==RealType(T)(0) || !A.isherm());
    TMVAssert(x.colsize() > 0);
    TMVAssert(x.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);

    if (x.isrm()) {
      if (A.isrm()) RowRankKUpdate(alpha,x,A);
      else ColRankKUpdate(alpha,x,A);
    } else {
      for (size_t i=0;i<x.rowsize();++i)
	Rank1Update(alpha,x.col(i),A);
    }
  }

#ifdef BLAS
  template <class T, class Tx> inline void BlasRankKUpdate(
      const T alpha, const GenMatrix<Tx>& x, const SymMatrixView<T>& A)
  { NonBlasRankKUpdate(alpha,x,A); }
  template <> inline void BlasRankKUpdate(
      const double alpha, const GenMatrix<double>& x,
      const SymMatrixView<double>& A)
  {
#ifdef XDEBUG
    //cerr<<"Blas RankK double\n";
#endif
    TMVAssert(A.size() == x.colsize());
    TMVAssert(alpha != double(0));
    TMVAssert(x.colsize() > 0);
    TMVAssert(x.rowsize() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.sym() == Sym);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(x.isrm() || x.iscm());

    if (A.isrm())
      if (x.isrm())
	cblas_dsyrk(CblasRowMajor,CblasLower,CblasNoTrans,
	    A.size(),x.rowsize(),alpha,
	    x.cptr(),x.stepi(),double(1),A.ptr(),A.stepi()); 
      else
	cblas_dsyrk(CblasRowMajor,CblasLower,CblasTrans,
	    A.size(),x.rowsize(),alpha,
	    x.cptr(),x.stepj(),double(1),A.ptr(),A.stepi()); 
    else
      if (x.isrm())
	cblas_dsyrk(CblasColMajor,CblasLower,CblasTrans,
	    A.size(),x.rowsize(),alpha,
	    x.cptr(),x.stepi(),double(1),A.ptr(),A.stepj()); 
      else
	cblas_dsyrk(CblasColMajor,CblasLower,CblasNoTrans,
	    A.size(),x.rowsize(),alpha,
	    x.cptr(),x.stepj(),double(1),A.ptr(),A.stepj()); 
  }
  template <> void BlasRankKUpdate(
      const complex<double> alpha, const GenMatrix<complex<double> >& x, 
      const SymMatrixView<complex<double> >& A)
  {
#ifdef XDEBUG
    //cerr<<"Blas RankK c double\n";
#endif
    TMVAssert(A.size() == x.colsize());
    TMVAssert(alpha != double(0));
    TMVAssert(IMAG(alpha)==double(0) || !A.isherm());
    TMVAssert(x.colsize() > 0);
    TMVAssert(x.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(x.isrm() || x.iscm());

    if (A.isherm()) {
      TMVAssert(IMAG(alpha)==double(0));
      if (A.isrm())
	if (x.isrm())
	  if (x.isconj()) NonBlasRankKUpdate(alpha,x,A);
	  else 
	    cblas_zherk(CblasRowMajor,CblasLower,CblasNoTrans,
		A.size(),x.rowsize(),REAL(alpha),
		x.cptr(),x.stepi(),double(1),A.ptr(),A.stepi()); 
	else
	  if (x.isconj())
	    cblas_zherk(CblasRowMajor,CblasLower,CblasConjTrans,
		A.size(),x.rowsize(),REAL(alpha),
		x.cptr(),x.stepj(),double(1),A.ptr(),A.stepi()); 
	  else NonBlasRankKUpdate(alpha,x,A);
      else
	if (x.isrm())
	  if (x.isconj())
	    cblas_zherk(CblasColMajor,CblasLower,CblasConjTrans,
		A.size(),x.rowsize(),REAL(alpha),
		x.cptr(),x.stepi(),double(1),A.ptr(),A.stepj()); 
	  else NonBlasRankKUpdate(alpha,x,A);
	else
	  if (x.isconj()) NonBlasRankKUpdate(alpha,x,A);
	  else 
	    cblas_zherk(CblasColMajor,CblasLower,CblasNoTrans,
		A.size(),x.rowsize(),REAL(alpha),
		x.cptr(),x.stepj(),double(1),A.ptr(),A.stepj()); 
    } else {
      if (x.isconj()) NonBlasRankKUpdate(alpha,x,A);
      else {
	complex<double> beta(1);
	if (A.isrm())
	  if (x.isrm())
	    cblas_zsyrk(CblasRowMajor,CblasLower,CblasNoTrans,
		A.size(),x.rowsize(),&alpha,
		x.cptr(),x.stepi(),&beta,A.ptr(),A.stepi()); 
	  else
	    cblas_zsyrk(CblasRowMajor,CblasLower,CblasTrans,
		A.size(),x.rowsize(),&alpha,
		x.cptr(),x.stepj(),&beta,A.ptr(),A.stepi()); 
	else
	  if (x.isrm())
	    cblas_zsyrk(CblasColMajor,CblasLower,CblasTrans,
		A.size(),x.rowsize(),&alpha,
		x.cptr(),x.stepi(),&beta,A.ptr(),A.stepj()); 
	  else
	    cblas_zsyrk(CblasColMajor,CblasLower,CblasNoTrans,
		A.size(),x.rowsize(),&alpha,
		x.cptr(),x.stepj(),&beta,A.ptr(),A.stepj()); 
      }
    }
  }
#ifndef NOFLOAT
  template <> inline void BlasRankKUpdate(
      const float alpha, const GenMatrix<float>& x,
      const SymMatrixView<float>& A)
  {
#ifdef XDEBUG
    //cerr<<"Blas RankK float\n";
#endif
    TMVAssert(A.size() == x.colsize());
    TMVAssert(alpha != float(0));
    TMVAssert(x.colsize() > 0);
    TMVAssert(x.rowsize() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.sym() == Sym);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(x.isrm() || x.iscm());

    if (A.isrm())
      if (x.isrm())
	cblas_ssyrk(CblasRowMajor,CblasLower,CblasNoTrans,
	    A.size(),x.rowsize(),alpha,
	    x.cptr(),x.stepi(),float(1),A.ptr(),A.stepi()); 
      else
	cblas_ssyrk(CblasRowMajor,CblasLower,CblasTrans,
	    A.size(),x.rowsize(),alpha,
	    x.cptr(),x.stepj(),float(1),A.ptr(),A.stepi()); 
    else
      if (x.isrm())
	cblas_ssyrk(CblasColMajor,CblasLower,CblasTrans,
	    A.size(),x.rowsize(),alpha,
	    x.cptr(),x.stepi(),float(1),A.ptr(),A.stepj()); 
      else
	cblas_ssyrk(CblasColMajor,CblasLower,CblasNoTrans,
	    A.size(),x.rowsize(),alpha,
	    x.cptr(),x.stepj(),float(1),A.ptr(),A.stepj()); 
  }
  template <> void BlasRankKUpdate(
      const complex<float> alpha, const GenMatrix<complex<float> >& x, 
      const SymMatrixView<complex<float> >& A)
  {
#ifdef XDEBUG
    //cerr<<"Blas RankK c float\n";
#endif
    TMVAssert(A.size() == x.colsize());
    TMVAssert(alpha != float(0));
    TMVAssert(IMAG(alpha)==float(0) || !A.isherm());
    TMVAssert(x.colsize() > 0);
    TMVAssert(x.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(x.isrm() || x.iscm());

    if (A.isherm()) {
      TMVAssert(IMAG(alpha)==float(0));
      if (A.isrm())
	if (x.isrm())
	  if (x.isconj()) NonBlasRankKUpdate(alpha,x,A);
	  else 
	    cblas_cherk(CblasRowMajor,CblasLower,CblasNoTrans,
		A.size(),x.rowsize(),REAL(alpha),
		x.cptr(),x.stepi(),float(1),A.ptr(),A.stepi()); 
	else
	  if (x.isconj())
	    cblas_cherk(CblasRowMajor,CblasLower,CblasConjTrans,
		A.size(),x.rowsize(),REAL(alpha),
		x.cptr(),x.stepj(),float(1),A.ptr(),A.stepi()); 
	  else NonBlasRankKUpdate(alpha,x,A);
      else
	if (x.isrm())
	  if (x.isconj())
	    cblas_cherk(CblasColMajor,CblasLower,CblasConjTrans,
		A.size(),x.rowsize(),REAL(alpha),
		x.cptr(),x.stepi(),float(1),A.ptr(),A.stepj()); 
	  else NonBlasRankKUpdate(alpha,x,A);
	else
	  if (x.isconj()) NonBlasRankKUpdate(alpha,x,A);
	  else 
	    cblas_cherk(CblasColMajor,CblasLower,CblasNoTrans,
		A.size(),x.rowsize(),REAL(alpha),
		x.cptr(),x.stepj(),float(1),A.ptr(),A.stepj()); 
    } else {
      if (x.isconj()) NonBlasRankKUpdate(alpha,x,A);
      else {
	complex<float> beta(1);
	if (A.isrm())
	  if (x.isrm())
	    cblas_csyrk(CblasRowMajor,CblasLower,CblasNoTrans,
		A.size(),x.rowsize(),&alpha,
		x.cptr(),x.stepi(),&beta,A.ptr(),A.stepi()); 
	  else
	    cblas_csyrk(CblasRowMajor,CblasLower,CblasTrans,
		A.size(),x.rowsize(),&alpha,
		x.cptr(),x.stepj(),&beta,A.ptr(),A.stepi()); 
	else
	  if (x.isrm())
	    cblas_csyrk(CblasColMajor,CblasLower,CblasTrans,
		A.size(),x.rowsize(),&alpha,
		x.cptr(),x.stepi(),&beta,A.ptr(),A.stepj()); 
	  else
	    cblas_csyrk(CblasColMajor,CblasLower,CblasNoTrans,
		A.size(),x.rowsize(),&alpha,
		x.cptr(),x.stepj(),&beta,A.ptr(),A.stepj()); 
      }
    }
  }
#endif // NOFLOAT
#endif // BLAS

  template <class T, class Tx> void RankKUpdate(
      const T alpha, const GenMatrix<Tx>& x, const SymMatrixView<T>& A)
    // A = A + alpha * x * xT
  {
#ifdef XDEBUG
    Matrix<T> A2 = Matrix<T>(A);
    if (A.isherm())
      A2 += (alpha*x^x.Conjugate());
    else 
      A2 += (alpha*x^x);
    Matrix<T> A0 = A;
    cerr<<"Start RankKUpdate: A = "<<Type(A)<<"  "<<A<<endl;
    cerr<<"alpha = "<<alpha<<", x = "<<Type(x)<<"  "<<x<<endl;
#endif

    TMVAssert(A.size() == x.colsize());
    TMVAssert(IMAG(alpha)==RealType(T)(0) || !A.isherm());
    if (alpha != T(0) && x.colsize() > 0 && x.rowsize() > 0) {
      if (A.isconj()) 
	RankKUpdate(CONJ(alpha),x.Conjugate(),A.Conjugate());
      else if (A.uplo() == Upper) {
	if (A.isherm()) RankKUpdate(alpha,x,A.Adjoint());
	else RankKUpdate(alpha,x,A.Transpose());
      }
#ifdef BLAS
      else if ((A.isrm() || A.iscm()) && (x.isrm() || x.iscm()))
	BlasRankKUpdate(alpha,x,A);
#endif
      else NonBlasRankKUpdate(alpha,x,A);
    }
  
#ifdef XDEBUG
    if (Norm(A2-A) > 0.001*Norm(A)) {
      cerr<<"RankKUpdate: alpha = "<<alpha<<endl;
      cerr<<"x = "<<Type(x)<<"  "<<x<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> A = "<<A<<endl;
      cerr<<"A2 = "<<A2<<endl;
      abort();
    }
#endif
  }

  // 
  // Rank2KUpdate
  //

  template <class T, class Tx, class Ty> void RowRank2KUpdate(
      const T alpha, const GenMatrix<Tx>& x,
      const GenMatrix<Ty>& y, const SymMatrixView<T>& A)
  {
#ifdef XDEBUG
    //cerr<<"Row Rank2K\n";
#endif
    TMVAssert(x.isrm());
    TMVAssert(A.isrm());
    TMVAssert(A.size() == x.colsize());
    TMVAssert(A.size() == y.colsize());
    TMVAssert(x.rowsize() == y.rowsize());
    TMVAssert(alpha != T(0));
    TMVAssert(x.colsize() > 0);
    TMVAssert(x.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);

    for (size_t i=0;i<A.size();++i) {
      if (A.isherm())
	A.row(i,0,i+1) += alpha * x.row(i) * y.Rows(0,i+1).Adjoint();
      else
	A.row(i,0,i+1) += alpha * x.row(i) * y.Rows(0,i+1).Transpose();
    }
  }

  template <class T, class Tx, class Ty> void ColRank2KUpdate(
      const T alpha, const GenMatrix<Tx>& x,
      const GenMatrix<Ty>& y, const SymMatrixView<T>& A)
  {
#ifdef XDEBUG
    //cerr<<"Col RankK\n";
#endif
    TMVAssert(y.isrm());
    TMVAssert(A.iscm());
    TMVAssert(A.size() == x.colsize());
    TMVAssert(A.size() == y.colsize());
    TMVAssert(x.rowsize() == y.rowsize());
    TMVAssert(alpha != T(0));
    TMVAssert(x.colsize() > 0);
    TMVAssert(x.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);

    const size_t N = A.size();

    for (size_t j=0;j<N;++j) {
      if (A.isherm())
	A.col(j,j,N) += alpha * x.Rows(j,N) * y.row(j).Conjugate();
      else
	A.col(j,j,N) += alpha * x.Rows(j,N) * x.row(j);
    }
  }

  template <class T, class Tx, class Ty> void NonBlasRank2KUpdate(
      const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
      const SymMatrixView<T>& A)
  {
    TMVAssert(A.size() == x.colsize());
    TMVAssert(A.size() == y.colsize());
    TMVAssert(x.rowsize() == y.rowsize());
    TMVAssert(alpha != T(0));
    TMVAssert(x.colsize() > 0);
    TMVAssert(x.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);

    if (x.isrm() && A.isrm()) RowRank2KUpdate(alpha,x,y,A);
    else if (y.isrm() && A.iscm()) ColRank2KUpdate(alpha,x,y,A);
    else
      for (size_t i=0;i<x.rowsize();i++)
	Rank2Update(alpha,x.col(i),y.col(i),A);
  }

#ifdef BLAS
  template <class T, class Tx, class Ty> inline void BlasRank2KUpdate(
      const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
      const SymMatrixView<T>& A)
  { NonBlasRank2KUpdate(alpha,x,y,A); }
  template <> inline void BlasRank2KUpdate(
      const double alpha, const GenMatrix<double>& x,
      const GenMatrix<double>& y, const SymMatrixView<double>& A)
  {
#ifdef XDEBUG
    //cerr<<"Blas Rank2K double\n";
#endif
    TMVAssert(A.size() == x.colsize());
    TMVAssert(A.size() == y.colsize());
    TMVAssert(x.rowsize() == y.rowsize());
    TMVAssert(alpha != double(0));
    TMVAssert(x.colsize() > 0);
    TMVAssert(x.rowsize() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.sym() == Sym);
    TMVAssert(x.isrm() || x.iscm());
    TMVAssert(x.stor() == y.stor());

    if (A.isrm())
      if (x.isrm())
	cblas_dsyr2k(CblasRowMajor,CblasLower,CblasNoTrans,
	    A.size(),x.rowsize(),alpha,
	    x.cptr(),x.stepi(),y.cptr(),y.stepi(),
	    double(1),A.ptr(),A.stepi()); 
      else
	cblas_dsyr2k(CblasRowMajor,CblasLower,CblasTrans,
	    A.size(),x.rowsize(),alpha,
	    x.cptr(),x.stepj(),y.cptr(),y.stepj(),
	    double(1),A.ptr(),A.stepi()); 
    else
      if (x.isrm())
	cblas_dsyr2k(CblasColMajor,CblasLower,CblasTrans,
	    A.size(),x.rowsize(),alpha,
	    x.cptr(),x.stepi(),y.cptr(),y.stepi(),
	    double(1),A.ptr(),A.stepj()); 
      else
	cblas_dsyr2k(CblasColMajor,CblasLower,CblasNoTrans,
	    A.size(),x.rowsize(),alpha,
	    x.cptr(),x.stepj(),y.cptr(),y.stepj(),
	    double(1),A.ptr(),A.stepj()); 
  }
  template <> void BlasRank2KUpdate(
      const complex<double> alpha, const GenMatrix<complex<double> >& x, 
      const GenMatrix<complex<double> >& y, 
      const SymMatrixView<complex<double> >& A)
  {
#ifdef XDEBUG
    //cerr<<"Blas Rank2K c double\n";
#endif
    TMVAssert(A.size() == x.colsize());
    TMVAssert(A.size() == y.colsize());
    TMVAssert(x.rowsize() == y.rowsize());
    TMVAssert(alpha != double(0));
    TMVAssert(x.colsize() > 0);
    TMVAssert(x.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.sym() == Sym);
    TMVAssert(x.isrm() || x.iscm());
    TMVAssert(x.stor() == y.stor());

    if (A.isherm()) {
      if (A.isrm())
	if (x.isrm()) 
	  if (x.isconj() || y.isconj()) NonBlasRank2KUpdate(alpha,x,y,A);
	  else 
	    cblas_zher2k(CblasRowMajor,CblasLower,CblasNoTrans,
		A.size(),x.rowsize(),&alpha,
		x.cptr(),x.stepi(),y.cptr(),y.stepi(),
		double(1),A.ptr(),A.stepi()); 
	else
	  if (x.isconj() && y.isconj()) 
	    cblas_zher2k(CblasRowMajor,CblasLower,CblasConjTrans,
		A.size(),x.rowsize(),&alpha,
		y.cptr(),y.stepj(),x.cptr(),x.stepj(),
		double(1),A.ptr(),A.stepi()); 
	  else 
	    NonBlasRank2KUpdate(alpha,x,y,A);
      else
	if (x.isrm())
	  if (x.isconj() && y.isconj()) 
	    cblas_zher2k(CblasColMajor,CblasLower,CblasConjTrans,
		A.size(),x.rowsize(),&alpha,
		y.cptr(),y.stepi(),x.cptr(),x.stepi(),
		double(1),A.ptr(),A.stepj()); 
	  else 
	    NonBlasRank2KUpdate(alpha,x,y,A);
	else
	  if (x.isconj() || y.isconj()) NonBlasRank2KUpdate(alpha,x,y,A);
	  else 
	    cblas_zher2k(CblasColMajor,CblasLower,CblasNoTrans,
		A.size(),x.rowsize(),&alpha,
		x.cptr(),x.stepj(),y.cptr(),y.stepj(),
		double(1),A.ptr(),A.stepj()); 
    }
    else {
      if (x.isconj() || y.isconj()) NonBlasRank2KUpdate(alpha,x,y,A);
      else {
	complex<double> beta(1);
	if (A.isrm())
	  if (x.isrm())
	    cblas_zsyr2k(CblasRowMajor,CblasLower,CblasNoTrans,
		A.size(),x.rowsize(),&alpha,
		x.cptr(),x.stepi(),y.cptr(),y.stepi(),
		&beta,A.ptr(),A.stepi()); 
	  else
	    cblas_zsyr2k(CblasRowMajor,CblasLower,CblasTrans,
		A.size(),x.rowsize(),&alpha,
		x.cptr(),x.stepj(),y.cptr(),y.stepj(),
		&beta,A.ptr(),A.stepi()); 
	else
	  if (x.isrm())
	    cblas_zsyr2k(CblasColMajor,CblasLower,CblasTrans,
		A.size(),x.rowsize(),&alpha,
		x.cptr(),x.stepi(),y.cptr(),y.stepi(),
		&beta,A.ptr(),A.stepj()); 
	  else
	    cblas_zsyr2k(CblasColMajor,CblasLower,CblasNoTrans,
		A.size(),x.rowsize(),&alpha,
		x.cptr(),x.stepj(),y.cptr(),y.stepj(),
		&beta,A.ptr(),A.stepj()); 
      }
    }
  }
#ifndef NOFLOAT
  template <> inline void BlasRank2KUpdate(
      const float alpha, const GenMatrix<float>& x,
      const GenMatrix<float>& y, const SymMatrixView<float>& A)
  {
#ifdef XDEBUG
    //cerr<<"Blas Rank2K float\n";
#endif
    TMVAssert(A.size() == x.colsize());
    TMVAssert(A.size() == y.colsize());
    TMVAssert(x.rowsize() == y.rowsize());
    TMVAssert(alpha != float(0));
    TMVAssert(x.colsize() > 0);
    TMVAssert(x.rowsize() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(y.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.sym() == Sym);
    TMVAssert(x.isrm() || x.iscm());
    TMVAssert(x.stor() == y.stor());

    if (A.isrm())
      if (x.isrm())
	cblas_ssyr2k(CblasRowMajor,CblasLower,CblasNoTrans,
	    A.size(),x.rowsize(),alpha,
	    x.cptr(),x.stepi(),y.cptr(),y.stepi(),
	    float(1),A.ptr(),A.stepi()); 
      else
	cblas_ssyr2k(CblasRowMajor,CblasLower,CblasTrans,
	    A.size(),x.rowsize(),alpha,
	    x.cptr(),x.stepj(),y.cptr(),y.stepj(),
	    float(1),A.ptr(),A.stepi()); 
    else
      if (x.isrm())
	cblas_ssyr2k(CblasColMajor,CblasLower,CblasTrans,
	    A.size(),x.rowsize(),alpha,
	    x.cptr(),x.stepi(),y.cptr(),y.stepi(),
	    float(1),A.ptr(),A.stepj()); 
      else
	cblas_ssyr2k(CblasColMajor,CblasLower,CblasNoTrans,
	    A.size(),x.rowsize(),alpha,
	    x.cptr(),x.stepj(),y.cptr(),y.stepj(),
	    float(1),A.ptr(),A.stepj()); 
  }
  template <> void BlasRank2KUpdate(
      const complex<float> alpha, const GenMatrix<complex<float> >& x, 
      const GenMatrix<complex<float> >& y, 
      const SymMatrixView<complex<float> >& A)
  {
#ifdef XDEBUG
    //cerr<<"Blas Rank2K c float\n";
#endif
    TMVAssert(A.size() == x.colsize());
    TMVAssert(A.size() == y.colsize());
    TMVAssert(x.rowsize() == y.rowsize());
    TMVAssert(alpha != float(0));
    TMVAssert(x.colsize() > 0);
    TMVAssert(x.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.sym() == Sym);
    TMVAssert(x.isrm() || x.iscm());
    TMVAssert(x.stor() == y.stor());

    if (A.isherm()) {
      if (A.isrm())
	if (x.isrm()) 
	  if (x.isconj() || y.isconj()) NonBlasRank2KUpdate(alpha,x,y,A);
	  else 
	    cblas_cher2k(CblasRowMajor,CblasLower,CblasNoTrans,
		A.size(),x.rowsize(),&alpha,
		x.cptr(),x.stepi(),y.cptr(),y.stepi(),
		float(1),A.ptr(),A.stepi()); 
	else
	  if (x.isconj() && y.isconj()) 
	    cblas_cher2k(CblasRowMajor,CblasLower,CblasConjTrans,
		A.size(),x.rowsize(),&alpha,
		y.cptr(),y.stepj(),x.cptr(),x.stepj(),
		float(1),A.ptr(),A.stepi()); 
	  else 
	    NonBlasRank2KUpdate(alpha,x,y,A);
      else
	if (x.isrm())
	  if (x.isconj() && y.isconj()) 
	    cblas_cher2k(CblasColMajor,CblasLower,CblasConjTrans,
		A.size(),x.rowsize(),&alpha,
		y.cptr(),y.stepi(),x.cptr(),x.stepi(),
		float(1),A.ptr(),A.stepj()); 
	  else 
	    NonBlasRank2KUpdate(alpha,x,y,A);
	else
	  if (x.isconj() || y.isconj()) NonBlasRank2KUpdate(alpha,x,y,A);
	  else 
	    cblas_cher2k(CblasColMajor,CblasLower,CblasNoTrans,
		A.size(),x.rowsize(),&alpha,
		x.cptr(),x.stepj(),y.cptr(),y.stepj(),
		float(1),A.ptr(),A.stepj()); 
    }
    else {
      if (x.isconj() || y.isconj()) NonBlasRank2KUpdate(alpha,x,y,A);
      else {
	complex<double> beta(1);
	if (A.isrm())
	  if (x.isrm())
	    cblas_csyr2k(CblasRowMajor,CblasLower,CblasNoTrans,
		A.size(),x.rowsize(),&alpha,
		x.cptr(),x.stepi(),y.cptr(),y.stepi(),
		&beta,A.ptr(),A.stepi()); 
	  else
	    cblas_csyr2k(CblasRowMajor,CblasLower,CblasTrans,
		A.size(),x.rowsize(),&alpha,
		x.cptr(),x.stepj(),y.cptr(),y.stepj(),
		&beta,A.ptr(),A.stepi()); 
	else
	  if (x.isrm())
	    cblas_csyr2k(CblasColMajor,CblasLower,CblasTrans,
		A.size(),x.rowsize(),&alpha,
		x.cptr(),x.stepi(),y.cptr(),y.stepi(),
		&beta,A.ptr(),A.stepj()); 
	  else
	    cblas_csyr2k(CblasColMajor,CblasLower,CblasNoTrans,
		A.size(),x.rowsize(),&alpha,
		x.cptr(),x.stepj(),y.cptr(),y.stepj(),
		&beta,A.ptr(),A.stepj()); 
      }
    }
  }
#endif // NOFLOAT
#endif // BLAS

  template <class T, class Tx, class Ty> void Rank2KUpdate(
      const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
      const SymMatrixView<T>& A)
    // if A is sym:  A = A + alpha * (x ^ y + y ^ x)
    // if A is herm: A = A + alpha * x ^ y* + conj(alpha) * y ^ x*
  {
#ifdef XDEBUG
    Matrix<T> A2 = Matrix<T>(A);
    if (A.isherm())
      A2 += (alpha*x^y.Conjugate())+(CONJ(alpha)*y^x.Conjugate());
    else
      A2 += alpha*((x^y) + (y^x));
    Matrix<T> A0 = A;
#endif

    TMVAssert(A.size() == x.colsize());
    TMVAssert(A.size() == y.rowsize());
    if (alpha != T(0) && A.size() > 0) {
      if (A.isconj()) 
	Rank2KUpdate(CONJ(alpha),x.Conjugate(),y.Conjugate(),
	    A.Conjugate());
      else if (A.uplo() == Upper) {
	if (A.isherm()) Rank2KUpdate(alpha,x,y,A.Adjoint());
	else Rank2KUpdate(alpha,x,y,A.Transpose());
      }
#ifdef BLAS
      else if ((A.isrm() || A.iscm()) && (x.isrm() || y.isrm()) &&
	  x.stor() == y.stor())
	BlasRank2KUpdate(alpha,x,y,A);
#endif
      else NonBlasRank2KUpdate(alpha,x,y,A);
    }
  
#ifdef XDEBUG
    if (Norm(A2-A) > 0.001*Norm(A)) {
      cerr<<"Rank2KUpdate: alpha = "<<alpha<<endl;
      cerr<<"x = "<<Type(x)<<"  "<<x<<endl;
      cerr<<"y = "<<Type(y)<<"  "<<y<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> A = "<<A<<endl;
      cerr<<"A2 = "<<A2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_SymMatrixArith_E.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


