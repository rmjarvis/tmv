
#include "TMV.h"
#include "TMV_Sym.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define SYM_RK_BLOCKSIZE TMV_BLOCKSIZE
#define SYM_RK_BLOCKSIZE2 (TMV_BLOCKSIZE/2)
#else
#define SYM_RK_BLOCKSIZE 64
#define SYM_RK_BLOCKSIZE2 1
#endif

  // 
  // RankKUpdate
  //

  template <bool ha, bool a1, bool add, class T, class Tx> 
    inline void RecursiveRankKUpdate(
	const T alpha, const GenMatrix<Tx>& x, const SymMatrixView<T>& A)
    {
      TMVAssert(A.size() == x.colsize());
      TMVAssert(alpha != T(0));
      TMVAssert(IMAG(alpha)==RealType(T)(0) || !A.isherm());
      TMVAssert(x.colsize() > 0);
      TMVAssert(x.rowsize() > 0);
      TMVAssert(A.ct() == NonConj);
      TMVAssert(A.uplo() == Lower);
      TMVAssert(ha == A.isherm());
      TMVAssert(a1 == (alpha == T(1)));

      size_t N = A.size();

      if (N <= SYM_RK_BLOCKSIZE2) {
	if (N == 1) {
	  if (ha) {
	    RealType(T) temp = x.row(0).NormSq();
	    TMVAssert(IMAG(alpha) == RealType(T)(0));
	    if (!a1) temp *= REAL(alpha);
	    if (add) *(A.ptr()) += temp;
	    else *(A.ptr()) = temp;
	  } else {
	    T temp = x.row(0) * x.row(0);
	    if (!a1) temp *= alpha;
	    if (add) *(A.ptr()) += temp;
	    else *(A.ptr()) = temp;
	  }
	} else {
	  if (x.isrm()) {
	    if (A.isrm()) {
	      for (size_t i=0;i<N;++i) {
		if (add)
		  A.row(i,0,i+1) += alpha * x.row(i) * 
		    (ha ? x.Rows(0,i+1).Adjoint() : x.Rows(0,i+1).Transpose());
		else
		  A.row(i,0,i+1) = alpha * x.row(i) * 
		    (ha ? x.Rows(0,i+1).Adjoint() : x.Rows(0,i+1).Transpose());
	      }
	    } else {
	      for (size_t j=0;j<N;++j) {
		if (add)
		  A.col(j,j,N) += alpha * x.Rows(j,N) * 
		    (ha ? x.row(j).Conjugate() : x.row(j));
		else
		  A.col(j,j,N) = alpha * x.Rows(j,N) * 
		    (ha ? x.row(j).Conjugate() : x.row(j));
	      }
	    }
	  } else { // x not row major
	    for (size_t i=0;i<x.rowsize();i++)
	      Rank1Update(alpha,x.col(i),add?1:0,A);
	  }
	}
      } else {
	size_t k = N/2;
	const size_t nb = SYM_RK_BLOCKSIZE;
	if (k > nb) k = k/nb*nb;
	RecursiveRankKUpdate<ha,a1,add>(alpha,x.Rows(0,k),A.SubSymMatrix(0,k));
	MultMM(alpha,x.Rows(k,N),
	    (ha ? x.Rows(0,k).Adjoint() : x.Rows(0,k).Transpose()),
	    add ? T(1) : T(0), A.SubMatrix(k,N,0,k));
	RecursiveRankKUpdate<ha,a1,add>(alpha,x.Rows(k,N),A.SubSymMatrix(k,N));
      }
    }

  template <bool cm, bool ha, bool a1, bool add, class T, class Tx> 
    inline void RecursiveInPlaceRankKUpdate(
	const T alpha, const GenMatrix<Tx>& x, const SymMatrixView<T>& A)
    {
      TMVAssert(x.Real().cptr() == A.Real().cptr());
      TMVAssert(A.size() > 0);
      TMVAssert(A.uplo() == Lower);
      TMVAssert(A.ct() == NonConj);
      TMVAssert(cm == A.iscm());
      TMVAssert(ha == A.isherm());
      TMVAssert(a1 == (alpha == T(1)));

      size_t N = A.size();
      if (N == 1) {
	Tx x00 = x(0,0);
	if (ha) {
	  RealType(T) temp = NORM(x00);
	  TMVAssert(IMAG(alpha) == RealType(T)(0));
	  if (!a1) temp *= REAL(alpha);
	  if (add) *A.ptr() += temp;
	  else *A.ptr() = temp;
	} else {
	  T temp = x00 * x00;
	  if (!a1) temp *= alpha;
	  if (add) *A.ptr() += temp;
	  else *A.ptr() = temp;
	}
      } else {
	// [ A00  A01 ] += alpha * [ x00  x01 ] [ x00t  x10t ]
	// [ A10  A11 ]            [ x10  x11 ] [ x01t  x11t ]
	//               = alpha * [ x00 x00t + x01 x01t   x00 x10t + x01 x11t ]
	//                         [ x10 x00t + x11 x01t   x10 x10t + x11 x11t ]
	// Note that there is no order to do these in which overwriting A??
	// won't screw up an x?? which is needed later.
	// Need a temporary.  I choose A10, but it doesn't much matter.
	const size_t k = N/2;
	const ConstMatrixView<Tx> x00 = x.SubMatrix(0,k,0,k);
	const ConstMatrixView<Tx> x10 = x.SubMatrix(k,N,0,k);
	const ConstMatrixView<Tx> x01 = x.SubMatrix(0,k,k,N);
	const ConstMatrixView<Tx> x11 = x.SubMatrix(k,N,k,N);
	SymMatrixView<T> A00 = A.SubSymMatrix(0,k);
	SymMatrixView<T> A11 = A.SubSymMatrix(k,N);
	MatrixView<T> A10 = A.SubMatrix(k,N,0,k);

	Matrix<T,cm?ColMajor:RowMajor> tempA10 = 
	  x10 * (ha ? x00.Adjoint() : x00.Transpose());
	tempA10 += x11 * (ha ? x01.Adjoint() : x01.Transpose());

	RecursiveInPlaceRankKUpdate<cm,ha,a1,add>(alpha,x11,A11);
	RecursiveRankKUpdate<ha,a1,true>(alpha,x10,A11);
	RecursiveInPlaceRankKUpdate<cm,ha,a1,add>(alpha,x00,A00);
	RecursiveRankKUpdate<ha,a1,true>(alpha,x01,A00);

	if (add) A10 += alpha * tempA10;
	else A10 = alpha * tempA10;
      }
    }

  template <bool add, class T, class Tx> inline void InPlaceRankKUpdate(
      const T alpha, const GenMatrix<Tx>& x, const SymMatrixView<T>& A)
  {
    TMVAssert(A.uplo() == Lower);
    if (A.iscm()) 
      if (A.isherm())
	if (alpha == T(1))
	  RecursiveInPlaceRankKUpdate<true,true,true,add>(alpha,x,A); 
	else
	  RecursiveInPlaceRankKUpdate<true,true,false,add>(alpha,x,A); 
      else
	if (alpha == T(1))
	  RecursiveInPlaceRankKUpdate<true,false,true,add>(alpha,x,A); 
	else
	  RecursiveInPlaceRankKUpdate<true,false,false,add>(alpha,x,A); 
    else
      if (A.isherm())
	if (alpha == T(1))
	  RecursiveInPlaceRankKUpdate<false,true,true,add>(alpha,x,A); 
	else
	  RecursiveInPlaceRankKUpdate<false,true,false,add>(alpha,x,A); 
      else
	if (alpha == T(1))
	  RecursiveInPlaceRankKUpdate<false,false,true,add>(alpha,x,A); 
	else
	  RecursiveInPlaceRankKUpdate<false,false,false,add>(alpha,x,A); 
  }

  template <bool add, class T, class Tx> inline void NonBlasRankKUpdate(
      const T alpha, const GenMatrix<Tx>& x, const SymMatrixView<T>& A)
  {
    if (A.uplo() == Upper) {
      return NonBlasRankKUpdate<add>(alpha,x,A.issym()?A.Transpose():A.Adjoint());
    } else if (A.isconj()) {
      return NonBlasRankKUpdate<add>(CONJ(alpha),x.Conjugate(),A.Conjugate());
    } else if (x.Real().cptr() == A.Real().cptr()) {
      const size_t N = A.size();
      TMVAssert(x.colsize() == N);
      const size_t K = x.rowsize();
      if (K >= N) {
	InPlaceRankKUpdate<add>(alpha,x.Cols(0,N),A);
	if (K > N) NonBlasRankKUpdate<add>(alpha,x.Cols(N,K),A);
      } else {
        NonBlasRankKUpdate<add>(alpha,x.Rows(K,N),A.SubSymMatrix(K,N));
	MultMM(alpha,x.Rows(K,N),
	    A.isherm() ? x.Rows(0,K).Adjoint() : x.Rows(0,K).Transpose(),
	    add ? T(1) : T(0), A.SubMatrix(K,N,0,K) );
	InPlaceRankKUpdate<add>(alpha,x.Rows(0,K),A.SubSymMatrix(0,K));
      }
    } else {
      if (A.isherm())
	if (alpha == T(1))
	  RecursiveRankKUpdate<true,true,add>(alpha,x,A); 
	else
	  RecursiveRankKUpdate<true,false,add>(alpha,x,A); 
      else
	if (alpha == T(1))
	  RecursiveRankKUpdate<false,true,add>(alpha,x,A); 
	else
	  RecursiveRankKUpdate<false,false,add>(alpha,x,A); 
    }
  }

#ifdef BLAS
  template <class T, class Tx> inline void BlasRankKUpdate(
      const T alpha, const GenMatrix<Tx>& x, int beta,
      const SymMatrixView<T>& A)
  { 
    if (beta == 0) NonBlasRankKUpdate<false>(alpha,x,A); 
    else NonBlasRankKUpdate<true>(alpha,x,A); 
  }
  template <> inline void BlasRankKUpdate(
      const double alpha, const GenMatrix<double>& x,
      int beta, const SymMatrixView<double>& A)
  {
    TMVAssert(A.size() == x.colsize());
    TMVAssert(alpha != double(0));
    TMVAssert(x.colsize() > 0);
    TMVAssert(x.rowsize() > 1);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.sym() == Sym);
    TMVAssert(A.iscm());
    TMVAssert(x.isrm() || x.iscm());

    int n = A.size();
    int k = x.rowsize();
    int ldx = x.isrm()?x.stepi():x.stepj();
    double b(beta);
    int lda = A.stepj();
    BLASNAME(dsyrk) (BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
	x.isrm()?BLASCH_T:BLASCH_NT,BLASV(n),BLASV(k),BLASV(alpha),
	BLASP(x.cptr()),BLASV(ldx),BLASV(b),BLASP(A.ptr()),BLASV(lda)
	BLAS1 BLAS1);
  }
  template <> inline void BlasRankKUpdate(
      const complex<double> alpha, const GenMatrix<complex<double> >& x, 
      int beta, const SymMatrixView<complex<double> >& A)
  {
    TMVAssert(A.size() == x.colsize());
    TMVAssert(alpha != double(0));
    TMVAssert(IMAG(alpha)==double(0) || !A.isherm());
    TMVAssert(x.colsize() > 0);
    TMVAssert(x.rowsize() > 1);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.iscm());
    TMVAssert(x.isrm() || x.iscm());
    TMVAssert(A.isherm() || !x.isconj());
    TMVAssert(A.issym() || x.isrm() == x.isconj());

    int n = A.size();
    int k = x.rowsize();
    int ldx = x.isrm()?x.stepi():x.stepj();
    int lda = A.stepj();
    if (A.isherm()) {
      double a(REAL(alpha));
      double b(beta);
      BLASNAME(zherk) (BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
	  x.isrm()?BLASCH_CT:BLASCH_NT,BLASV(n),BLASV(k),BLASV(a),
	  BLASP(x.cptr()),BLASV(ldx),BLASV(b),BLASP(A.ptr()),BLASV(lda)
	  BLAS1 BLAS1);
    } else {
      complex<double> b(beta);
      BLASNAME(zsyrk) (BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
	  x.isrm()?BLASCH_T:BLASCH_NT,BLASV(n),BLASV(k),BLASP(&alpha),
	  BLASP(x.cptr()),BLASV(ldx),BLASP(&b),BLASP(A.ptr()),BLASV(lda)
	  BLAS1 BLAS1);
    }
  }
#ifndef NOFLOAT
  template <> inline void BlasRankKUpdate(
      const float alpha, const GenMatrix<float>& x,
      int beta, const SymMatrixView<float>& A)
  {
    TMVAssert(A.size() == x.colsize());
    TMVAssert(alpha != float(0));
    TMVAssert(x.colsize() > 0);
    TMVAssert(x.rowsize() > 1);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.sym() == Sym);
    TMVAssert(A.iscm());
    TMVAssert(x.isrm() || x.iscm());

    int n = A.size();
    int k = x.rowsize();
    int ldx = x.isrm()?x.stepi():x.stepj();
    float b(beta);
    int lda = A.stepj();
    BLASNAME(ssyrk) (BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
	x.isrm()?BLASCH_T:BLASCH_NT,BLASV(n),BLASV(k),BLASV(alpha),
	BLASP(x.cptr()),BLASV(ldx),BLASV(b),BLASP(A.ptr()),BLASV(lda)
	BLAS1 BLAS1);
  }
  template <> inline void BlasRankKUpdate(
      const complex<float> alpha, const GenMatrix<complex<float> >& x, 
      int beta, const SymMatrixView<complex<float> >& A)
  {
    TMVAssert(A.size() == x.colsize());
    TMVAssert(alpha != float(0));
    TMVAssert(IMAG(alpha)==float(0) || !A.isherm());
    TMVAssert(x.colsize() > 0);
    TMVAssert(x.rowsize() > 1);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.iscm());
    TMVAssert(x.isrm() || x.iscm());
    TMVAssert(A.isherm() || !x.isconj());
    TMVAssert(A.issym() || x.isrm() == x.isconj());

    int n = A.size();
    int k = x.rowsize();
    int ldx = x.isrm()?x.stepi():x.stepj();
    int lda = A.stepj();
    if (A.isherm()) {
      float a(REAL(alpha));
      float b(beta);
      BLASNAME(cherk) (BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
	  x.isrm()?BLASCH_CT:BLASCH_NT,BLASV(n),BLASV(k),BLASV(a),
	  BLASP(x.cptr()),BLASV(ldx),BLASV(b),BLASP(A.ptr()),BLASV(lda)
	  BLAS1 BLAS1);
    } else {
      complex<float> b(beta);
      BLASNAME(csyrk) (BLASCM A.uplo()==Upper?BLASCH_UP:BLASCH_LO,
	  x.isrm()?BLASCH_T:BLASCH_NT,BLASV(n),BLASV(k),BLASP(&alpha),
	  BLASP(x.cptr()),BLASV(ldx),BLASP(&b),BLASP(A.ptr()),BLASV(lda)
	  BLAS1 BLAS1);
    }
  }
#endif // NOFLOAT
#endif // BLAS

  template <class T, class Tx> void RankKUpdate(
      const T alpha, const GenMatrix<Tx>& x, const int beta,
      const SymMatrixView<T>& A)
    // A = beta*A + alpha * x * xT
  {
#ifdef XDEBUG
    Matrix<T> A0 = A;
    Matrix<T> x0 = x;
    Matrix<T> A2 = Matrix<T>(A);
    Matrix<T> xxt = Matrix<T>(x)*Matrix<T>(
	A.isherm() ? x.Adjoint() : x.Transpose());
    if (beta == 0) A2 = alpha*xxt;
    else A2 += alpha * xxt;
    //cerr<<"Start RankKUpdate: A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"alpha = "<<alpha<<", x = "<<Type(x)<<"  "<<x<<endl;
#endif

    TMVAssert(A.size() == x.colsize());
    TMVAssert(beta == 0 || beta == 1);
    TMVAssert(IMAG(alpha)==RealType(T)(0) || !A.isherm());
    if (alpha != T(0) && x.colsize() > 0 && x.rowsize() > 0) {
      if (x.rowsize() == 1)
	return Rank1Update(alpha,x.col(0),beta,A);
#ifdef BLAS
      else if (IsComplex(T()) && IsReal(Tx()))
	BlasRankKUpdate(alpha,x,beta,A);
      else if (A.isrm())
	return RankKUpdate(alpha,x,beta,A.issym()?A.Transpose():A.Adjoint());
      else if (A.isconj()) 
	return RankKUpdate(CONJ(alpha),x.Conjugate(),beta,A.Conjugate());
      else if (A.iscm() && A.stepj()>0) {
	if (!((x.iscm() && x.stepj()>0) || (x.isrm() && x.stepi()>0)) || 
	    (!A.issym() && x.iscm() == x.isconj()) ||
	    (!A.isherm() && x.isconj()) ||
	    (x.Real().cptr() == A.Real().ptr())) {
	  Matrix<T,ColMajor> xx = x;
	  BlasRankKUpdate(alpha,xx,beta,A);
	} else {
	  BlasRankKUpdate(alpha,x,beta,A);
	}
      } else {
	if (A.isherm()) {
	  if (beta == 0) {
	    HermMatrix<T,Lower,ColMajor> AA(A.size(),RealType(T)(0));
	    RankKUpdate(alpha,x,beta,AA.View());
	    A = AA;
	  } else {
	    HermMatrix<T,Lower,ColMajor> AA = A;
	    RankKUpdate(alpha,x,beta,AA.View());
	    A = AA;
	  } 
	} else {
	  if (beta == 0) {
	    SymMatrix<T,Lower,ColMajor> AA(A.size(),T(0));
	    RankKUpdate(alpha,x,beta,AA.View());
	    A = AA;
	  } else {
	    SymMatrix<T,Lower,ColMajor> AA = A;
	    RankKUpdate(alpha,x,beta,AA.View());
	    A = AA;
	  }
	}
      } 
#else
      else if (beta == 0) NonBlasRankKUpdate<false>(alpha,x,A);
      else NonBlasRankKUpdate<true>(alpha,x,A);
#endif
    }
  
#ifdef XDEBUG
    if (Norm(A-A2) > 0.001*(abs(alpha)*SQR(Norm(x0))+Norm(A0))) {
      cerr<<"RankKUpdate: alpha = "<<alpha<<endl;
      cerr<<"x = "<<Type(x)<<"  "<<x0<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> A = "<<A<<endl;
      cerr<<"A2 = "<<A2<<endl;
      abort();
    }
#endif
  }

  // 
  // S += U*Ut
  //

  template <bool ha, bool uu, bool a1, class T, class TU> inline void RecursiveRankKUpdate(
      const T alpha, const GenUpperTriMatrix<TU>& U, const SymMatrixView<T>& A)
  {
    TMVAssert(A.size() == U.size());
    TMVAssert(alpha != T(0));
    TMVAssert(IMAG(alpha)==RealType(T)(0) || !A.isherm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Upper);
    TMVAssert(ha == A.isherm());
    TMVAssert(uu == U.isunit());
    TMVAssert(a1 == (alpha == T(1)));

    size_t N = A.size();

    if (N == 1) {
      if (uu)
	*(A.ptr()) += a1 ? T(1) : alpha;
      else 
	if (a1)
	  *(A.ptr()) += ha ? NORM(*(U.cptr())) : SQR(*(U.cptr()));
	else
	  *(A.ptr()) += alpha * (ha ? NORM(*(U.cptr())) : SQR(*(U.cptr())));
    } else {
      //
      // [ A11 A12 ] += alpha [ U11 U12 ] [ U11t  0   ]
      // [ A21 A22 ]          [  0  U22 ] [ U12t U22t ]
      //              = alpha [ U11 U11t + U12 U12t    U12 U22t ]
      //                      [       U22 U12t         U22 U22t ]
      size_t k = N/2;
      const size_t nb = SYM_RK_BLOCKSIZE;
      if (k > nb) k = k/nb*nb;
      SymMatrixView<T> A11 = A.SubSymMatrix(0,k);
      SymMatrixView<T> A22 = A.SubSymMatrix(k,N);
      MatrixView<T> A12 = A.SubMatrix(0,k,k,N);
      ConstUpperTriMatrixView<TU> U11 = U.SubTriMatrix(0,k);
      ConstUpperTriMatrixView<TU> U22 = U.SubTriMatrix(k,N);
      ConstMatrixView<TU> U12 = U.SubMatrix(0,k,k,N);

      RecursiveRankKUpdate<ha,uu,a1>(alpha,U11,A11);

      RankKUpdate(alpha,U12,1,A11);

      A12 += alpha * U12 * (ha ? U22.Adjoint() : U22.Transpose());

      RecursiveRankKUpdate<ha,uu,a1>(alpha,U22,A22);
    }
  }

  template <bool a1, class T, class TU> inline void DoRankKUpdate(
      const T alpha, const GenUpperTriMatrix<TU>& U, const SymMatrixView<T>& A)
  {
    if (A.isherm())
      if (U.isunit())
	RecursiveRankKUpdate<true,true,a1>(alpha,U,A);
      else
	RecursiveRankKUpdate<true,false,a1>(alpha,U,A);
    else
      if (U.isunit())
	RecursiveRankKUpdate<false,true,a1>(alpha,U,A);
      else
	RecursiveRankKUpdate<false,false,a1>(alpha,U,A);
  }

  template <bool ha, class T> inline void RecursiveSetUUt(const SymMatrixView<T>& A)
  {
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Upper);
    TMVAssert(ha == A.isherm());

    size_t N = A.size();

    if (N == 1) {
      *(A.ptr()) = ha ? NORM(*(A.cptr())) : SQR(*(A.cptr()));
    } else {
      //
      // [ A11 A12 ] = alpha [ U11 U12 ] [ U11t  0   ]
      // [ A21 A22 ]         [  0  U22 ] [ U12t U22t ]
      //             = alpha [ U11 U11t + U12 U12t    U12 U22t ]
      //                     [       U22 U12t         U22 U22t ]
      size_t k = N/2;
      const size_t nb = SYM_RK_BLOCKSIZE;
      if (k > nb) k = k/nb*nb;
      SymMatrixView<T> A11 = A.SubSymMatrix(0,k);
      SymMatrixView<T> A22 = A.SubSymMatrix(k,N);
      MatrixView<T> A12 = A.SubMatrix(0,k,k,N);
      ConstUpperTriMatrixView<T> U22 = A22.UpperTri();

      RecursiveSetUUt<ha>(A11);

      // The Transpose's here are because these RankKUpdate routines
      // want A11 to be stored in the Lower triangle.
      if (ha)
	RankKUpdate(T(1),A12.Conjugate(),1,A11.Transpose());
      else
	RankKUpdate(T(1),A12,1,A11.Transpose());

      A12 *= (ha ? U22.Adjoint() : U22.Transpose());

      RecursiveSetUUt<ha>(A22);
    }
  }

#ifdef AELAP
  template <class T> inline void LapSetUUt(const SymMatrixView<T>& A)
  { RecursiveSetUUt<false>(A); }
  template<> inline void LapSetUUt(const SymMatrixView<double>& A)
  {
    TMVAssert(A.issym());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.uplo() == Upper);

    int n = A.size();
    int lda = A.isrm() ? A.stepi() : A.stepj();
    LAPNAME(dlauum) (LAPCM A.isrm()?LAPCH_LO:LAPCH_UP,LAPV(n),
	LAPP(A.ptr()),LAPV(lda) LAPINFO LAP1);
    LAP_Results("dlauum");
  }
  template<> inline void LapSetUUt(const SymMatrixView<complex<double> >& A)
  {
    TMVAssert(A.issym());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.uplo() == Upper);

    int n = A.size();
    int lda = A.isrm() ? A.stepi() : A.stepj();
    LAPNAME(zlauum) (LAPCM A.isrm()?LAPCH_LO:LAPCH_UP,LAPV(n),
	LAPP(A.ptr()),LAPV(lda) LAPINFO LAP1);
    LAP_Results("zlauum");
  }
#ifndef NOFLOAT
  template<> inline void LapSetUUt(const SymMatrixView<float>& A)
  {
    TMVAssert(A.issym());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.uplo() == Upper);

    int n = A.size();
    int lda = A.isrm() ? A.stepi() : A.stepj();
    LAPNAME(slauum) (LAPCM A.isrm()?LAPCH_LO:LAPCH_UP,LAPV(n),
	LAPP(A.ptr()),LAPV(lda) LAPINFO LAP1);
    LAP_Results("slauum");
  }
  template<> inline void LapSetUUt(const SymMatrixView<complex<float> >& A)
  {
    TMVAssert(A.issym());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.uplo() == Upper);

    int n = A.size();
    int lda = A.isrm() ? A.stepi() : A.stepj();
    LAPNAME(clauum) (LAPCM A.isrm()?LAPCH_LO:LAPCH_UP,LAPV(n),
	LAPP(A.ptr()),LAPV(lda) LAPINFO LAP1);
    LAP_Results("clauum");
  }
#endif // FLOAT
#endif // AELAP

  template <class T> inline void SetUUt(const SymMatrixView<T>& A)
  {
#ifdef AELAP
    if (A.issym() && ((A.isrm() && A.stepi()>0) || (A.iscm() && A.stepj()>0)))
      LapSetUUt(A);
    else
#endif
      if (A.isherm()) RecursiveSetUUt<true>(A);
      else RecursiveSetUUt<false>(A);
  }

  template <class T, class TU> void RankKUpdate(
      const T alpha, const GenUpperTriMatrix<TU>& U, const int beta,
      const SymMatrixView<T>& A)
    // A = A + alpha * U * UT
  {
#ifdef XDEBUG
    Matrix<T> A0 = A;
    Matrix<TU> U0 = U;
    Matrix<T> A2 = Matrix<T>(A);
    if (beta == 0)
      A2 = alpha*U*(A.isherm() ? U.Adjoint() : U.Transpose());
    else 
      A2 += alpha*U*(A.isherm() ? U.Adjoint() : U.Transpose());
    //cerr<<"Start RankKUpdate: A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"alpha = "<<alpha<<", U = "<<Type(U)<<"  "<<U<<endl;
#endif

    TMVAssert(A.size() == U.size());
    TMVAssert(beta == 0 || beta == 1);
    TMVAssert(IMAG(alpha)==RealType(T)(0) || !A.isherm());

    if (alpha != T(0) && A.size() > 0) {
      if (A.isconj()) 
	RankKUpdate(CONJ(alpha),U.Conjugate(),beta,A.Conjugate());
      else if (A.uplo() == Lower) 
	if (A.issym()) RankKUpdate(alpha,U,beta,A.Transpose());
	else RankKUpdate(CONJ(alpha),U.Conjugate(),beta,A.Transpose());
      else if (beta == 0) {
	A.UpperTri() = U;
	SetUUt(A);
	if (alpha != T(1)) A *= alpha;
      }
      else
	if (alpha == T(1)) DoRankKUpdate<true>(alpha,U,A);
	else DoRankKUpdate<false>(alpha,U,A);
    }
  
#ifdef XDEBUG
    if (Norm(A-A2) > 0.001*(abs(alpha)*SQR(Norm(U0))+Norm(A0))) {
      cerr<<"RankKUpdate: alpha = "<<alpha<<endl;
      cerr<<"U = "<<Type(U)<<"  "<<U0<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> A = "<<A<<endl;
      cerr<<"A2 = "<<A2<<endl;
      abort();
    }
#endif
  }

  // 
  // S += L*Lt
  //

  template <bool ha, bool uu, bool a1, class T, class TL> inline void RecursiveRankKUpdate(
      const T alpha, const GenLowerTriMatrix<TL>& L, const SymMatrixView<T>& A)
  {
    TMVAssert(A.size() == L.size());
    TMVAssert(alpha != T(0));
    TMVAssert(IMAG(alpha)==RealType(T)(0) || !A.isherm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(ha == A.isherm());
    TMVAssert(uu == L.isunit());
    TMVAssert(a1 == (alpha == T(1)));

    size_t N = A.size();

    if (N == 1) {
      if (uu)
	*(A.ptr()) += a1 ? T(1) : alpha;
      else 
	if (a1)
	  *(A.ptr()) += ha ? NORM(*(L.cptr())) : SQR(*(L.cptr()));
	else
	  *(A.ptr()) += alpha * (ha ? NORM(*(L.cptr())) : SQR(*(L.cptr())));
    } else {
      //
      // [ A11 A12 ] += alpha [ L11  0  ] [ L11t L21t ]
      // [ A21 A22 ]          [ L21 L22 ] [  0   L22t ]
      //              = alpha [ L11 L11t        L11 L21t      ]
      //                      [ L21 L11t  L21 L21t + L22 L22t ]
      size_t k = N/2;
      const size_t nb = SYM_RK_BLOCKSIZE;
      if (k > nb) k = k/nb*nb;
      SymMatrixView<T> A11 = A.SubSymMatrix(0,k);
      SymMatrixView<T> A22 = A.SubSymMatrix(k,N);
      MatrixView<T> A21 = A.SubMatrix(k,N,0,k);
      ConstLowerTriMatrixView<TL> L11 = L.SubTriMatrix(0,k);
      ConstLowerTriMatrixView<TL> L22 = L.SubTriMatrix(k,N);
      ConstMatrixView<TL> L21 = L.SubMatrix(k,N,0,k);

      RecursiveRankKUpdate<ha,uu,a1>(alpha,L22,A22);

      RankKUpdate(alpha,L21,1,A22);

      A21 += alpha * L21 * (ha ? L11.Adjoint() : L11.Transpose());

      RecursiveRankKUpdate<ha,uu,a1>(alpha,L11,A11);
    }
  }

  template <bool a1, class T, class TL> inline void DoRankKUpdate(
      const T alpha, const GenLowerTriMatrix<TL>& L, const SymMatrixView<T>& A)
  {
    if (A.isherm())
      if (L.isunit())
	RecursiveRankKUpdate<true,true,a1>(alpha,L,A);
      else
	RecursiveRankKUpdate<true,false,a1>(alpha,L,A);
    else
      if (L.isunit())
	RecursiveRankKUpdate<false,true,a1>(alpha,L,A);
      else
	RecursiveRankKUpdate<false,false,a1>(alpha,L,A);
  }

  template <bool ha, class T> inline void RecursiveSetLLt(const SymMatrixView<T>& A)
  {
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.uplo() == Lower);
    TMVAssert(ha == A.isherm());

    size_t N = A.size();

    if (N == 1) {
      *(A.ptr()) = ha ? NORM(*(A.cptr())) : SQR(*(A.cptr()));
    } else {
      // [ A11 A12 ] = [ L11  0  ] [ L11t L21t ]
      // [ A21 A22 ]   [ L21 L22 ] [  0   L22t ]
      //             = [ L11 L11t        L11 L21t      ]
      //               [ L21 L11t  L21 L21t + L22 L22t ]
      size_t k = N/2;
      const size_t nb = SYM_RK_BLOCKSIZE;
      if (k > nb) k = k/nb*nb;
      SymMatrixView<T> A11 = A.SubSymMatrix(0,k);
      SymMatrixView<T> A22 = A.SubSymMatrix(k,N);
      MatrixView<T> A21 = A.SubMatrix(k,N,0,k);
      ConstLowerTriMatrixView<T> L11 = A11.LowerTri();

      RecursiveSetLLt<ha>(A22);

      RankKUpdate(T(1),A21,1,A22);

      A21 *= (ha ? L11.Adjoint() : L11.Transpose());

      RecursiveSetLLt<ha>(A11);
    }
  }

  template <class T> inline void SetLLt(const SymMatrixView<T>& A)
  {
    if (A.isherm()) RecursiveSetLLt<true>(A);
    else RecursiveSetLLt<false>(A);
  }

  template <class T, class TL> void RankKUpdate(
      const T alpha, const GenLowerTriMatrix<TL>& L, const int beta,
      const SymMatrixView<T>& A)
    // A = A + alpha * L * LT
  {
#ifdef XDEBUG
    Matrix<T> A0 = A;
    Matrix<TL> L0 = L;
    Matrix<T> A2 = Matrix<T>(A);
    if (beta == 0)
      A2 = alpha*L*(A.isherm() ? L.Adjoint() : L.Transpose());
    else 
      A2 += alpha*L*(A.isherm() ? L.Adjoint() : L.Transpose());
    //cerr<<"Start RankKUpdate: A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"alpha = "<<alpha<<", L = "<<Type(L)<<"  "<<L<<endl;
#endif

    TMVAssert(A.size() == L.size());
    TMVAssert(beta == 0 || beta == 1);
    TMVAssert(IMAG(alpha)==RealType(T)(0) || !A.isherm());

    if (alpha != T(0) && A.size() > 0) {
      if (A.isconj()) 
	RankKUpdate(CONJ(alpha),L.Conjugate(),beta,A.Conjugate());
      else if (A.uplo() == Upper) 
	if (A.isherm()) RankKUpdate(alpha,L,beta,A.Adjoint());
	else RankKUpdate(alpha,L,beta,A.Transpose());
      else if (beta == 0) {
	A.LowerTri() = L;
	SetLLt(A);
	if (alpha != T(1)) A *= alpha;
      }
      else
	if (alpha == T(1)) DoRankKUpdate<true>(alpha,L,A);
	else DoRankKUpdate<false>(alpha,L,A);
    }
  
#ifdef XDEBUG
    if (Norm(A-A2) > 0.001*(abs(alpha)*SQR(Norm(L0))+
	  (beta==0?RealType(T)(0):Norm(A0)))) {
      cerr<<"RankKUpdate: alpha = "<<alpha<<endl;
      cerr<<"L = "<<Type(L)<<"  "<<L0<<endl;
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


