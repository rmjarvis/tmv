
#include "TMV.h"
#include "TMV_Sym.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
  const size_t SYM_R2K_BLOCKSIZE = TMV_BLOCKSIZE;
  const size_t SYM_R2K_BLOCKSIZE2 = TMV_BLOCKSIZE/2;
#else
  const size_t SYM_R2K_BLOCKSIZE = 64;
  const size_t SYM_R2K_BLOCKSIZE2 = 32;
#endif

  // 
  // Rank2KUpdate
  //

  template <class T, class Tx, class Ty> void NonBlasRank2KUpdate(
      const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
      const SymMatrixView<T>& A)
  {
    TMVAssert(A.size() == x.colsize());
    TMVAssert(A.size() == y.colsize());
    TMVAssert(x.rowsize() == y.rowsize());
    TMVAssert(alpha != T(0));
    TMVAssert(x.colsize() > 0);
    TMVAssert(x.rowsize() > 1);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);

    const size_t nb = SYM_R2K_BLOCKSIZE;
    size_t N = A.size();

    if (N <= SYM_R2K_BLOCKSIZE2) {
      if (x.isrm() && y.isrm()) {
	if (A.isrm()) {
	  for (size_t i=0;i<A.size();++i) {
	    if (A.isherm()) {
	      A.row(i,0,i+1) += alpha * x.row(i) * y.Rows(0,i+1).Adjoint();
	      A.row(i,0,i+1) += CONJ(alpha) * y.row(i) * 
		x.Rows(0,i+1).Adjoint();
	    } else {
	      A.row(i,0,i+1) += alpha * x.row(i) * y.Rows(0,i+1).Transpose();
	      A.row(i,0,i+1) += alpha * y.row(i) * x.Rows(0,i+1).Transpose();
	    }
	  }
	}
	else {
	  for (size_t j=0;j<N;++j) {
	    if (A.isherm())  {
	      A.col(j,j,N) += alpha * x.Rows(j,N) * y.row(j).Conjugate();
	      A.col(j,j,N) += CONJ(alpha) * y.Rows(j,N) * 
		x.row(j).Conjugate();
	    } else {
	      A.col(j,j,N) += alpha * x.Rows(j,N) * y.row(j);
	      A.col(j,j,N) += alpha * y.Rows(j,N) * x.row(j);
	    }
	  }
	}
      } else { // x,y not row major
	for (size_t i=0;i<x.rowsize();i++)
	  Rank2Update(alpha,x.col(i),y.col(i),A);
      }
    } else {
      size_t k = N/2;
      if (k > nb) k = k/nb*nb;
      NonBlasRank2KUpdate(alpha,x.Rows(0,k),y.Rows(0,k),A.SubSymMatrix(0,k));
      if (A.isherm()) {
	A.SubMatrix(k,N,0,k) += alpha * x.Rows(k,N) * y.Rows(0,k).Adjoint();
	A.SubMatrix(k,N,0,k) += CONJ(alpha) * y.Rows(k,N) * 
	  x.Rows(0,k).Adjoint();
      } else {
	A.SubMatrix(k,N,0,k) += alpha * x.Rows(k,N) * y.Rows(0,k).Transpose();
	A.SubMatrix(k,N,0,k) += alpha * y.Rows(k,N) * x.Rows(0,k).Transpose();
      }
      NonBlasRank2KUpdate(alpha,x.Rows(k,N),y.Rows(k,N),A.SubSymMatrix(k,N));
    }
  }

#ifdef BLAS
  template <class T, class Tx, class Ty> void BlasRank2KUpdate(
      const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
      const SymMatrixView<T>& A)
  { NonBlasRank2KUpdate(alpha,x,y,A); }
  template <> void BlasRank2KUpdate(
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
    TMVAssert(x.rowsize() > 1);
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
    TMVAssert(x.rowsize() > 1);
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
  template <> void BlasRank2KUpdate(
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
    TMVAssert(x.rowsize() > 1);
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
    TMVAssert(x.rowsize() > 1);
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
    TMVAssert(A.size() == x.colsize());
    TMVAssert(A.size() == y.colsize());
    TMVAssert(x.rowsize() == y.rowsize());

#ifdef XDEBUG
    Matrix<T> A2 = Matrix<T>(A);
    if (A.isherm())
      A2 += alpha*x*y.Adjoint()+CONJ(alpha)*y*x.Adjoint();
    else
      A2 += alpha*(x*y.Transpose() + y*x.Transpose());
    Matrix<T> A0 = A;
#endif

    if (alpha != T(0) && A.size() > 0) {
      if (A.isconj()) 
	Rank2KUpdate(CONJ(alpha),x.Conjugate(),y.Conjugate(),
	    A.Conjugate());
      else if (A.uplo() == Upper) {
	if (A.isherm()) Rank2KUpdate(alpha,x,y,A.Adjoint());
	else Rank2KUpdate(alpha,x,y,A.Transpose());
      }
      else if (x.rowsize() == 1)
	Rank2Update(alpha,x.col(0),y.col(0),A);
#ifdef BLAS
      else if ( ((A.isrm() && A.stepi()>0) || (A.iscm() && A.stepj()>0)) &&
	  ( (x.isrm() && y.isrm() && x.stepi()>0 && y.stepi()>0) || 
	    (x.iscm() && y.iscm() && x.stepj()>0 && y.stepj()>0) ) )
	BlasRank2KUpdate(alpha,x,y,A);
#endif
      else NonBlasRank2KUpdate(alpha,x,y,A);
    }
  
#ifdef XDEBUG
    if (Norm(A-A2) > 0.001*max(RealType(T)(1),Norm(A))) {
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

#define InstFile "TMV_SymMatrixArith_H.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


