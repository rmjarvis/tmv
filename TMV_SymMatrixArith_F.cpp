
#include "TMV.h"
#include "TMV_Sym.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define SYM_R2K_BLOCKSIZE TMV_BLOCKSIZE
#define SYM_R2K_BLOCKSIZE2 (TMV_BLOCKSIZE/2)
#else
#define SYM_R2K_BLOCKSIZE 64
#define SYM_R2K_BLOCKSIZE2 1
#endif

  // 
  // Rank2KUpdate
  //

  template <bool ha, bool a1, bool add, class T, class Tx, class Ty> 
    inline void RecursiveRank2KUpdate(
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
      TMVAssert(ha == A.isherm());
      TMVAssert(a1 == (alpha == T(1)));

      const size_t nb = SYM_R2K_BLOCKSIZE;
      size_t N = A.size();

      if (N <= SYM_R2K_BLOCKSIZE2) {
	if (N == 1) {
	  T temp = x.row(0) * (ha ? y.row(0).Conjugate() : y.row(0));
	  if (!a1) temp *= alpha;
	  if (ha)
	    if (add) *(A.ptr()) += RealType(T)(2) * REAL(temp);
	    else *(A.ptr()) = RealType(T)(2) * REAL(temp);
	  else
	    if (add) *(A.ptr()) += RealType(T)(2) * temp;
	    else *(A.ptr()) = RealType(T)(2) * temp;
	} else {
	  if (x.isrm() && y.isrm()) {
	    if (A.isrm()) {
	      for (size_t i=0;i<N;++i) {
		if (add) 
		  A.row(i,0,i+1) += alpha * x.row(i) * 
		    (ha ? y.Rows(0,i+1).Adjoint() : y.Rows(0,i+1).Transpose());
		else 
		  A.row(i,0,i+1) = alpha * x.row(i) * 
		    (ha ? y.Rows(0,i+1).Adjoint() : y.Rows(0,i+1).Transpose());
		A.row(i,0,i+1) += y.row(i) * 
		  (ha ? (CONJ(alpha)*x.Rows(0,i+1).Adjoint()) : 
		   (alpha*x.Rows(0,i+1).Transpose()));
	      }
	    }
	    else {
	      for (size_t j=0;j<N;++j) {
		if (add) 
		  A.col(j,j,N) += alpha * x.Rows(j,N) * 
		    (ha ? y.row(j).Conjugate() : y.row(j));
		else 
		  A.col(j,j,N) = alpha * x.Rows(j,N) * 
		    (ha ? y.row(j).Conjugate() : y.row(j));
		A.col(j,j,N) += y.Rows(j,N) * 
		  (ha ? (CONJ(alpha)*x.row(j).Conjugate()) : (alpha*x.row(j)));
	      }
	    }
	  } else { // x,y not row major
	    for (size_t i=0;i<x.rowsize();++i) {
	      Rank2Update(alpha,x.col(i),y.col(i),add?1:0,A);
	    }
	  }
	}
      } else { // Not <= BLOCKSIZE2, so do recurse...
	size_t k = N/2;
	if (k > nb) k = k/nb*nb;

	RecursiveRank2KUpdate<ha,a1,add>(alpha,x.Rows(0,k),y.Rows(0,k),
	    A.SubSymMatrix(0,k));

	if (add) 
	  A.SubMatrix(k,N,0,k) += alpha * x.Rows(k,N) * 
	    (ha ? y.Rows(0,k).Adjoint() : y.Rows(0,k).Transpose());
	else 
	  A.SubMatrix(k,N,0,k) = alpha * x.Rows(k,N) * 
	    (ha ? y.Rows(0,k).Adjoint() : y.Rows(0,k).Transpose());

	A.SubMatrix(k,N,0,k) += y.Rows(k,N) * 
	  (ha ? (CONJ(alpha)*x.Rows(0,k).Adjoint()) : 
	   (alpha * x.Rows(0,k).Transpose()));

	RecursiveRank2KUpdate<ha,a1,add>(alpha,x.Rows(k,N),y.Rows(k,N),
	    A.SubSymMatrix(k,N));
      }
    }

  template <bool ha, bool a1, bool add, class T, class Tx, class Ty> 
    inline void RecursiveInPlaceRank2KUpdate(
	const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
	const SymMatrixView<T>& A)
    {
      TMVAssert(x.Real().cptr() == A.Real().cptr() || 
	  y.Real().cptr() == A.Real().cptr());
      TMVAssert(A.size() > 0);
      TMVAssert(A.uplo() == Lower);
      TMVAssert(A.ct() == NonConj);

      size_t N = A.size();
      if (N == 1) {
	Tx x00 = x(0,0);
	Ty y00 = y(0,0);
	if (ha) {
	  RealType(T) temp = RealType(T)(2) *
	    (a1 ? REAL(x00*CONJ(y00)) : REAL(alpha*(x00*(CONJ(y00)))));
	  if (add) *A.ptr() += temp;
	  else *A.ptr() = temp;
	} else {
	  T temp = RealType(T)(2) * 
	    (a1 ? (x00 * y00) : (alpha * (x00 * y00)));
	  if (add) *A.ptr() += temp;
	  else *A.ptr() = temp;
	}
      } else {
	const size_t k = N/2;
	const ConstMatrixView<Tx> x00 = x.SubMatrix(0,k,0,k);
	const ConstMatrixView<Tx> x10 = x.SubMatrix(k,N,0,k);
	const ConstMatrixView<Tx> x01 = x.SubMatrix(0,k,k,N);
	const ConstMatrixView<Tx> x11 = x.SubMatrix(k,N,k,N);
	const ConstMatrixView<Ty> y00 = y.SubMatrix(0,k,0,k);
	const ConstMatrixView<Ty> y10 = y.SubMatrix(k,N,0,k);
	const ConstMatrixView<Ty> y01 = y.SubMatrix(0,k,k,N);
	const ConstMatrixView<Ty> y11 = y.SubMatrix(k,N,k,N);
	SymMatrixView<T> A00 = A.SubSymMatrix(0,k);
	SymMatrixView<T> A11 = A.SubSymMatrix(k,N);
	MatrixView<T> A10 = A.SubMatrix(k,N,0,k);

	Matrix<T> tempA10 = x10 * (ha ? y00.Adjoint() : y00.Transpose());
	tempA10 += x11 * (ha ? y01.Adjoint() : y01.Transpose());
	if (!a1) tempA10 *= alpha;
	T ca = ha ? CONJ(alpha) : alpha;
	tempA10 += ca * y10 * (ha ? x00.Adjoint() : x00.Transpose());
	tempA10 += ca * y11 * (ha ? x01.Adjoint() : x01.Transpose());
	RecursiveInPlaceRank2KUpdate<ha,a1,add>(alpha,x11,y11,A11);
	RecursiveRank2KUpdate<ha,a1,true>(alpha,x10,y10,A11);
	RecursiveInPlaceRank2KUpdate<ha,a1,add>(alpha,x00,y00,A00);
	RecursiveRank2KUpdate<ha,a1,true>(alpha,x01,y01,A00);

	if (add) A10 += tempA10;
	else A10 = tempA10;

      }
    }
							    
							    
  template <class T, class Tx, class Ty> inline void InPlaceRank2KUpdate(
      const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
      const int beta, const SymMatrixView<T>& A)
  {
    TMVAssert(x.Real().cptr() == A.Real().cptr() || 
	y.Real().cptr() == A.Real().cptr());
    TMVAssert(A.size() > 0);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.ct() == NonConj);
 
    if (A.isherm())
      if (alpha == T(1))
	if (beta == 0)
	  RecursiveInPlaceRank2KUpdate<true,true,false>(alpha,x,y,A); 
	else
	  RecursiveInPlaceRank2KUpdate<true,true,true>(alpha,x,y,A); 
      else
	if (beta == 0)
	  RecursiveInPlaceRank2KUpdate<true,false,false>(alpha,x,y,A); 
	else
	  RecursiveInPlaceRank2KUpdate<true,false,true>(alpha,x,y,A); 
    else
      if (alpha == T(1))
	if (beta == 0)
	  RecursiveInPlaceRank2KUpdate<false,true,false>(alpha,x,y,A); 
	else
	  RecursiveInPlaceRank2KUpdate<false,true,true>(alpha,x,y,A); 
      else
	if (beta == 0)
	  RecursiveInPlaceRank2KUpdate<false,false,false>(alpha,x,y,A); 
	else
	  RecursiveInPlaceRank2KUpdate<false,false,true>(alpha,x,y,A); 
  } 
  template <class T, class Tx, class Ty> inline void NonBlasRank2KUpdate(
      const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
      const int beta, const SymMatrixView<T>& A)
  { 
    if (A.uplo() == Upper)
      return NonBlasRank2KUpdate(alpha,x,y,beta,A.issym()?A.Transpose():A.Adjoint());
    else if (A.isconj())
      return NonBlasRank2KUpdate(CONJ(alpha),x.Conjugate(),y.Conjugate(),beta,
	  A.Conjugate());
    else if (x.Real().cptr() == A.Real().cptr() || 
	y.Real().cptr() == A.Real().cptr()) {
      const size_t N = A.size();
      TMVAssert(x.colsize() == N);
      TMVAssert(y.colsize() == N);
      const size_t K = x.rowsize();
      TMVAssert(y.rowsize() == K);

      if (K >= N) {
	InPlaceRank2KUpdate(alpha,x.Cols(0,N),y.Cols(0,N),beta,A);
	if (K > N) NonBlasRank2KUpdate(alpha,x.Cols(N,K),y.Cols(N,K),beta,A);
      } else { 
	NonBlasRank2KUpdate(alpha,x.Rows(K,N),y.Rows(K,N),beta,
	    A.SubSymMatrix(K,N));
	const bool ha = A.isherm();
	if (x.stepi() < x.stepj() == A.stepi() < A.stepj()) {
	  if (y.stepi() < y.stepj() == A.stepi() < A.stepj()) {
	    // Then both x and y overlap with A(K:N,0:K)
	    // Need a temporary
	    if (A.iscm()) {
	      if (ha) {
		Matrix<T,ColMajor> temp = 
		  alpha * x.Rows(K,N) * y.Rows(0,K).Adjoint();
		temp += CONJ(alpha) * y.Rows(K,N) * x.Rows(0,K).Adjoint();
		AddMM(T(1),temp,T(beta),A.SubMatrix(K,N,0,K));
	      } else {
		Matrix<T,ColMajor> temp = 
		  alpha * x.Rows(K,N) * y.Rows(0,K).Transpose();
		temp += alpha * y.Rows(K,N) * x.Rows(0,K).Transpose();
		AddMM(T(1),temp,T(beta),A.SubMatrix(K,N,0,K));
	      }
	    } else {
	      if (ha) {
		Matrix<T,RowMajor> temp = 
		  alpha * x.Rows(K,N) * y.Rows(0,K).Adjoint();
		temp += CONJ(alpha) * y.Rows(K,N) * x.Rows(0,K).Adjoint();
		AddMM(T(1),temp,T(beta),A.SubMatrix(K,N,0,K));
	      } else {
		Matrix<T,RowMajor> temp = 
		  alpha * x.Rows(K,N) * y.Rows(0,K).Transpose();
		temp += alpha * y.Rows(K,N) * x.Rows(0,K).Transpose();
		AddMM(T(1),temp,T(beta),A.SubMatrix(K,N,0,K));
	      }
	    }
	  } else {
	    MultMM(alpha, x.Rows(K,N),
		ha ? y.Rows(0,K).Adjoint() : y.Rows(0,K).Transpose(),
		T(1), A.SubMatrix(K,N,0,K) );
	    MultMM(ha ? CONJ(alpha) : alpha, y.Rows(K,N),
		ha ? x.Rows(0,K).Adjoint() : x.Rows(0,K).Transpose(),
		T(beta), A.SubMatrix(K,N,0,K) );
	  }
	}
	else {
	  MultMM(ha ? CONJ(alpha) : alpha, y.Rows(K,N),
	      ha ? x.Rows(0,K).Adjoint() : x.Rows(0,K).Transpose(),
	      T(1), A.SubMatrix(K,N,0,K) );
	  MultMM(alpha, x.Rows(K,N),
	      ha ? y.Rows(0,K).Adjoint() : y.Rows(0,K).Transpose(),
	      T(beta), A.SubMatrix(K,N,0,K) );
	}
	InPlaceRank2KUpdate(alpha,x.Rows(0,K),y.Rows(0,K),beta,
	    A.SubSymMatrix(0,K));
      }
    }
    else 
      if (A.isherm())
	if (alpha == T(1))
	  if (beta == 0)
	    RecursiveRank2KUpdate<true,true,false>(alpha,x,y,A); 
	  else
	    RecursiveRank2KUpdate<true,true,true>(alpha,x,y,A); 
	else
	  if (beta == 0)
	    RecursiveRank2KUpdate<true,false,false>(alpha,x,y,A); 
	  else
	    RecursiveRank2KUpdate<true,false,true>(alpha,x,y,A); 
      else
	if (alpha == T(1))
	  if (beta == 0)
	    RecursiveRank2KUpdate<false,true,false>(alpha,x,y,A); 
	  else
	    RecursiveRank2KUpdate<false,true,true>(alpha,x,y,A); 
	else
	  if (beta == 0)
	    RecursiveRank2KUpdate<false,false,false>(alpha,x,y,A); 
	  else
	    RecursiveRank2KUpdate<false,false,true>(alpha,x,y,A); 
  }

#ifdef BLAS
  template <class T, class Tx, class Ty> inline void BlasRank2KUpdate(
      const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
      const int beta, const SymMatrixView<T>& A)
  { NonBlasRank2KUpdate(alpha,x,y,beta,A); }
  template <> inline void BlasRank2KUpdate(
      const double alpha, const GenMatrix<double>& x,
      const GenMatrix<double>& y, const int beta,
      const SymMatrixView<double>& A)
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
    TMVAssert(A.iscm());
    TMVAssert(x.isrm() || x.iscm());
    TMVAssert(y.isrm() || y.iscm());
    TMVAssert(x.stor() == y.stor());

    int n=A.size();
    int k=x.rowsize();
    int ldx=x.isrm()?x.stepi():x.stepj();
    int ldy=y.isrm()?y.stepi():y.stepj();
    double b(beta);
    int lda=A.stepj();
    BLASNAME(dsyr2k) (BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO, 
	x.isrm()?BLASCH_T:BLASCH_NT,BLASV(n),BLASV(k),BLASV(alpha),
	BLASP(x.cptr()),BLASV(ldx),BLASP(y.cptr()),BLASV(ldy),
	BLASV(b),BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1);
  }
  template <> inline void BlasRank2KUpdate(
      const complex<double> alpha, const GenMatrix<complex<double> >& x, 
      const GenMatrix<complex<double> >& y, const int beta,
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
    TMVAssert(A.iscm());
    TMVAssert(x.isrm() || x.iscm());
    TMVAssert(y.isrm() || y.iscm());
    TMVAssert(x.stor() == y.stor());
    TMVAssert(x.isconj() == y.isconj());
    TMVAssert(A.isherm() || !x.isconj());
    TMVAssert(A.issym() || x.isrm() == x.isconj());

    int n=A.size();
    int k=x.rowsize();
    int ldx=x.isrm()?x.stepi():x.stepj();
    int ldy=y.isrm()?y.stepi():y.stepj();
    int lda=A.stepj();
    if (A.isherm()) {
      double b(beta);
      BLASNAME(zher2k) (BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO, 
	  x.isrm()?BLASCH_CT:BLASCH_NT,BLASV(n),BLASV(k),BLASP(&alpha),
	  BLASP(x.cptr()),BLASV(ldx),BLASP(y.cptr()),BLASV(ldy),
	  BLASV(b),BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1);
    }
    else {
      complex<double> b(beta);
      BLASNAME(zsyr2k) (BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO, 
	  x.isrm()?BLASCH_T:BLASCH_NT,BLASV(n),BLASV(k),BLASP(&alpha),
	  BLASP(x.cptr()),BLASV(ldx),BLASP(y.cptr()),BLASV(ldy),
	  BLASP(&b),BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1);
    }
  }
#ifndef NOFLOAT
  template <> inline void BlasRank2KUpdate(
      const float alpha, const GenMatrix<float>& x,
      const GenMatrix<float>& y, const int beta,
      const SymMatrixView<float>& A)
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
    TMVAssert(A.iscm());
    TMVAssert(x.isrm() || x.iscm());
    TMVAssert(y.isrm() || y.iscm());
    TMVAssert(x.stor() == y.stor());

    int n=A.size();
    int k=x.rowsize();
    int ldx=x.isrm()?x.stepi():x.stepj();
    int ldy=y.isrm()?y.stepi():y.stepj();
    float b(beta);
    int lda=A.stepj();
    BLASNAME(ssyr2k) (BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO, 
	x.isrm()?BLASCH_T:BLASCH_NT,BLASV(n),BLASV(k),BLASV(alpha),
	BLASP(x.cptr()),BLASV(ldx),BLASP(y.cptr()),BLASV(ldy),
	BLASV(b),BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1);
  }
  template <> inline void BlasRank2KUpdate(
      const complex<float> alpha, const GenMatrix<complex<float> >& x, 
      const GenMatrix<complex<float> >& y, const int beta,
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
    TMVAssert(A.iscm());
    TMVAssert(x.isrm() || x.iscm());
    TMVAssert(y.isrm() || y.iscm());
    TMVAssert(x.stor() == y.stor());
    TMVAssert(x.isconj() == y.isconj());
    TMVAssert(A.isherm() || !x.isconj());
    TMVAssert(A.issym() || x.isrm() == x.isconj());

    int n=A.size();
    int k=x.rowsize();
    int ldx=x.isrm()?x.stepi():x.stepj();
    int ldy=y.isrm()?y.stepi():y.stepj();
    int lda=A.stepj();
    if (A.isherm()) {
      float b(beta);
      BLASNAME(cher2k) (BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO, 
	  x.isrm()?BLASCH_CT:BLASCH_NT,BLASV(n),BLASV(k),BLASP(&alpha),
	  BLASP(x.cptr()),BLASV(ldx),BLASP(y.cptr()),BLASV(ldy),
	  BLASV(b),BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1);
    }
    else {
      complex<float> b(beta);
      BLASNAME(csyr2k) (BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO, 
	  x.isrm()?BLASCH_T:BLASCH_NT,BLASV(n),BLASV(k),BLASP(&alpha),
	  BLASP(x.cptr()),BLASV(ldx),BLASP(y.cptr()),BLASV(ldy),
	  BLASP(&b),BLASP(A.ptr()),BLASV(lda) BLAS1 BLAS1);
    }
  }
#endif // NOFLOAT
#endif // BLAS

  template <class T, class Tx, class Ty> void Rank2KUpdate(
      const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
      const int beta, const SymMatrixView<T>& A)
    // if A is sym:  A = beta*A + alpha * (x ^ y + y ^ x)
    // if A is herm: A = beta*A + alpha * x ^ y* + conj(alpha) * y ^ x*
  {
    TMVAssert(A.size() == x.colsize());
    TMVAssert(A.size() == y.colsize());
    TMVAssert(x.rowsize() == y.rowsize());

#ifdef XDEBUG
    Matrix<T> A0 = A;
    Matrix<Tx> x0 = x;
    Matrix<Ty> y0 = y;
    Matrix<T> A2(A);
    if (A.isherm()) {
      if (beta == 0) A2 = alpha*x*y.Adjoint();
      else A2 += alpha*x*y.Adjoint();
      A2 += CONJ(alpha)*y*x.Adjoint();
    }
    else {
      if (beta == 0) A2 = alpha*x*y.Transpose();
      else A2 += alpha*x*y.Transpose();
      A2 += alpha*y*x.Transpose();
    }
#endif

    TMVAssert(beta == 0 || beta == 1);
    if (alpha != T(0) && A.size() > 0) {
      if (x.rowsize() == 1)
	Rank2Update(alpha,x.col(0),y.col(0),beta,A);
#ifdef BLAS
      else if (IsComplex(T()) && (IsReal(Tx()) || IsReal(Ty())))
	BlasRank2KUpdate(alpha,x,y,beta,A);
      else if (A.isrm())
	Rank2KUpdate(alpha,x,y,beta,A.issym()?A.Transpose():A.Adjoint());
      else if (A.isconj()) 
	Rank2KUpdate(CONJ(alpha),x.Conjugate(),y.Conjugate(),beta,
	    A.Conjugate());
      else if (A.iscm() && A.stepj()>0) {
	if (!((x.isrm() && x.stepi()>0) || (x.iscm() && x.stepj()>0)) ||
	    (!A.issym() && x.iscm() == x.isconj()) ||
	    (!A.isherm() && x.isconj()) || 
	    x.Real().cptr() == A.Real().ptr()) {
	  if (!((y.isrm() && y.stepi()>0) || (y.iscm() && y.stepj()>0)) ||
	      y.isconj() ||
	      y.Real().cptr() == A.Real().ptr()) {
	    if (IMAG(alpha) == RealType(T)(0)) {
	      Matrix<Tx,ColMajor> xx = REAL(alpha)*x;
	      Matrix<Ty,ColMajor> yy = y;
	      BlasRank2KUpdate(T(1),xx,yy,beta,A);
	    } else {
	      Matrix<T,ColMajor> xx = alpha*x;
	      Matrix<Ty,ColMajor> yy = y;
	      BlasRank2KUpdate(T(1),xx,yy,beta,A);
	    }
	  } else {
	    if (y.iscm()) {
	      if (IMAG(alpha) == RealType(T)(0)) {
		Matrix<Tx,ColMajor> xx = REAL(alpha)*x;
		BlasRank2KUpdate(T(1),xx,y,beta,A);
	      } else {
		Matrix<T,ColMajor> xx = alpha*x;
		BlasRank2KUpdate(T(1),xx,y,beta,A);
	      }
	    } else {
	      if (IMAG(alpha) == RealType(T)(0)) {
		Matrix<Tx,RowMajor> xx = REAL(alpha)*x;
		BlasRank2KUpdate(T(1),xx,y,beta,A);
	      } else {
		Matrix<T,RowMajor> xx = alpha*x;
		BlasRank2KUpdate(T(1),xx,y,beta,A);
	      }
	    }
	  }
	} else {
	  if (!((y.isrm() && y.stepi()>0) || (y.iscm() && y.stepj()>0)) ||
	      y.stor() != x.stor() || y.isconj() ||
	      y.Real().cptr() == A.Real().ptr()) {
	    if (x.iscm()) {
	      if (IMAG(alpha) == RealType(T)(0)) {
		Matrix<Ty,ColMajor> yy = REAL(alpha)*y;
		BlasRank2KUpdate(T(1),x,yy,beta,A);
	      } else {
		Matrix<T,ColMajor> yy = alpha*y;
		BlasRank2KUpdate(T(1),x,yy,beta,A);
	      }
	    } else {
	      if (IMAG(alpha) == RealType(T)(0)) {
		Matrix<Ty,RowMajor> yy = REAL(alpha)*y;
		BlasRank2KUpdate(T(1),x,yy,beta,A);
	      } else {
		Matrix<T,RowMajor> yy = alpha*y;
		BlasRank2KUpdate(T(1),x,yy,beta,A);
	      }
	    }
	  } else {
	    BlasRank2KUpdate(alpha,x,y,beta,A);
	  }
	}
      } else {
	if (A.isherm()) {
	  if (beta == 0) {
	    HermMatrix<T,Lower,ColMajor> AA(A.size(),RealType(T)(0));
	    Rank2KUpdate(alpha,x,y,1,AA.View());
	    A = AA;
	  } else {
	    HermMatrix<T,Lower,ColMajor> AA = A;
	    Rank2KUpdate(alpha,x,y,1,AA.View());
	    A = AA;
	  }
	} else {
	  if (beta == 0) {
	    SymMatrix<T,Lower,ColMajor> AA(A.size(),T(0));
	    Rank2KUpdate(alpha,x,y,1,AA.View());
	    A = AA;
	  } else {
	    SymMatrix<T,Lower,ColMajor> AA = A;
	    Rank2KUpdate(alpha,x,y,1,AA.View());
	    A = AA;
	  }
	}
      }
#else
      else NonBlasRank2KUpdate(alpha,x,y,beta,A);
#endif
    }
  
#ifdef XDEBUG
    if (Norm(A-A2) > 0.001*(abs(alpha)*Norm(x0)*Norm(y0)+
	  (beta==0?RealType(T)(0):Norm(A0)))) {
      cerr<<"Rank2KUpdate: alpha,beta = "<<alpha<<','<<beta<<endl;
      cerr<<"x = "<<Type(x)<<"  "<<x0<<endl;
      cerr<<"y = "<<Type(y)<<"  "<<y0<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      if (A.isherm())
	cerr<<"x*yt = "<<x0*y0.Adjoint()<<endl;
      else
	cerr<<"x*yT = "<<x0*y0.Transpose()<<endl;
      cerr<<"-> A = "<<A<<endl;
      cerr<<"A2 = "<<A2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_SymMatrixArith_F.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


