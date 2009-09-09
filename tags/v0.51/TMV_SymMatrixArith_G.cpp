
#include "TMV.h"
#include "TMV_Sym.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
  const size_t SYM_RK_BLOCKSIZE = TMV_BLOCKSIZE;
  const size_t SYM_RK_BLOCKSIZE2 = TMV_BLOCKSIZE/2;
#else
  const size_t SYM_RK_BLOCKSIZE = 64;
  const size_t SYM_RK_BLOCKSIZE2 = 32;
#endif

  // 
  // RankKUpdate
  //

  template <class T, class Tx> void NonBlasRankKUpdate(
      const T alpha, const GenMatrix<Tx>& x, const SymMatrixView<T>& A)
  {
    TMVAssert(A.size() == x.colsize());
    TMVAssert(alpha != T(0));
    TMVAssert(IMAG(alpha)==RealType(T)(0) || !A.isherm());
    TMVAssert(x.colsize() > 0);
    TMVAssert(x.rowsize() > 1);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);

    const size_t nb = SYM_RK_BLOCKSIZE;
    size_t N = A.size();

    if (N <= SYM_RK_BLOCKSIZE2) {
      if (x.isrm()) {
	if (A.isrm()) {
	  for (size_t i=0;i<A.size();++i) {
	    if (A.isherm()) 
	      A.row(i,0,i+1) += alpha * x.row(i) * x.Rows(0,i+1).Adjoint();
	    else 
	      A.row(i,0,i+1) += alpha * x.row(i) * x.Rows(0,i+1).Transpose();
	  }
	} else {
	  for (size_t j=0;j<N;++j) {
	    if (A.isherm()) 
	      A.col(j,j,N) += alpha * x.Rows(j,N) * x.row(j).Conjugate();
	    else 
	      A.col(j,j,N) += alpha * x.Rows(j,N) * x.row(j);
	  }
	}
      } else { // x not row major
	for (size_t i=0;i<x.rowsize();i++)
	  Rank1Update(alpha,x.col(i),A);
      }
    } else {
      size_t k = N/2;
      if (k > nb) k = k/nb*nb;
      NonBlasRankKUpdate(alpha,x.Rows(0,k),A.SubSymMatrix(0,k));
      if (A.issym()) 
	A.SubMatrix(k,N,0,k) += alpha * x.Rows(k,N) * x.Rows(0,k).Transpose();
      else 
	A.SubMatrix(k,N,0,k) += alpha * x.Rows(k,N) * x.Rows(0,k).Adjoint();
      NonBlasRankKUpdate(alpha,x.Rows(k,N),A.SubSymMatrix(k,N));
    }
  }

#ifdef BLAS
  template <class T, class Tx> void BlasRankKUpdate(
      const T alpha, const GenMatrix<Tx>& x, const SymMatrixView<T>& A)
  { NonBlasRankKUpdate(alpha,x,A); }
  template <> void BlasRankKUpdate(
      const double alpha, const GenMatrix<double>& x,
      const SymMatrixView<double>& A)
  {
#ifdef XDEBUG
    //cerr<<"Blas RankK double\n";
#endif
    TMVAssert(A.size() == x.colsize());
    TMVAssert(alpha != double(0));
    TMVAssert(x.colsize() > 0);
    TMVAssert(x.rowsize() > 1);
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
    TMVAssert(x.rowsize() > 1);
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
  template <> void BlasRankKUpdate(
      const float alpha, const GenMatrix<float>& x,
      const SymMatrixView<float>& A)
  {
#ifdef XDEBUG
    //cerr<<"Blas RankK float\n";
#endif
    TMVAssert(A.size() == x.colsize());
    TMVAssert(alpha != float(0));
    TMVAssert(x.colsize() > 0);
    TMVAssert(x.rowsize() > 1);
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
    TMVAssert(x.rowsize() > 1);
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
      A2 += alpha*x*x.Adjoint();
    else 
      A2 += alpha*x*x.Transpose();
    Matrix<T> A0 = A;
    //cerr<<"Start RankKUpdate: A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"alpha = "<<alpha<<", x = "<<Type(x)<<"  "<<x<<endl;
#endif

    TMVAssert(A.size() == x.colsize());
    TMVAssert(IMAG(alpha)==RealType(T)(0) || !A.isherm());
    if (alpha != T(0) && x.colsize() > 0 && x.rowsize() > 0) {
      if (A.isconj()) 
	RankKUpdate(CONJ(alpha),x.Conjugate(),A.Conjugate());
      else if (A.uplo() == Upper) 
	if (A.isherm()) RankKUpdate(alpha,x,A.Adjoint());
	else RankKUpdate(alpha,x,A.Transpose());
      else if (x.rowsize() == 1)
	Rank1Update(alpha,x.col(0),A);
#ifdef BLAS
      else if ( ((A.isrm() && A.stepi()>0) || (A.iscm() && A.stepj()>0)) 
	  && ((x.isrm() && x.stepi()>0) || (x.iscm() && x.stepj()>0)))
	BlasRankKUpdate(alpha,x,A);
#endif
      else NonBlasRankKUpdate(alpha,x,A);
    }
  
#ifdef XDEBUG
    if (Norm(A-A2) > 0.001*max(RealType(T)(1),Norm(A))) {
      cerr<<"RankKUpdate: alpha = "<<alpha<<endl;
      cerr<<"x = "<<Type(x)<<"  "<<x<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> A = "<<A<<endl;
      cerr<<"A2 = "<<A2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_SymMatrixArith_G.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


