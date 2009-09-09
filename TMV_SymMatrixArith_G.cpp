
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
  // SymMultMM
  // A += alpha * x * y
  // where x,y are Matrices, and the product x*y is assumed to be symmetric.
  //

  template <bool ha, bool a1, bool add, class T, class Tx, class Ty> 
    inline void RecursiveSymMultMM(
	const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
	const SymMatrixView<T>& A)
    {
      TMVAssert(A.size() == x.colsize());
      TMVAssert(A.size() == y.rowsize());
      TMVAssert(x.rowsize() == y.colsize());
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
	  T temp = x.row(0) * y.col(0);
	  if (!a1) temp *= alpha;
	  if (ha)
	    if (add) *(A.ptr()) += REAL(temp);
	    else *(A.ptr()) = REAL(temp);
	  else
	    if (add) *(A.ptr()) += temp;
	    else *(A.ptr()) = temp;
	} else {
	  if (A.isrm()) {
	    for (size_t i=0;i<N;++i) {
	      if (add) 
		A.row(i,0,i+1) += alpha * x.row(i) * y.Cols(0,i+1);
	      else 
		A.row(i,0,i+1) = alpha * x.row(i) * y.Cols(0,i+1);
	    }
	  } else {
	    for (size_t j=0;j<N;++j) {
	      if (add) 
		A.col(j,j,N) += alpha * x.Rows(j,N) * y.col(j);
	      else 
		A.col(j,j,N) = alpha * x.Rows(j,N) * y.col(j);
	    }
	  }
	}
      } else { // Not <= BLOCKSIZE2, so do recurse...
	size_t k = N/2;
	if (k > nb) k = k/nb*nb;

	RecursiveSymMultMM<ha,a1,add>(alpha,x.Rows(0,k),y.Cols(0,k),
	    A.SubSymMatrix(0,k));

	if (add) 
	  A.SubMatrix(k,N,0,k) += alpha * x.Rows(k,N) * y.Cols(0,k);
	else 
	  A.SubMatrix(k,N,0,k) = alpha * x.Rows(k,N) * y.Cols(0,k);

	RecursiveSymMultMM<ha,a1,add>(alpha,x.Rows(k,N),y.Cols(k,N),
	    A.SubSymMatrix(k,N));
      }
    }

    template <bool ha, bool a1, bool add, class T, class Tx, class Ty> 
      inline void RecursiveInPlaceSymMultMM(
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
	    RealType(T) temp = a1 ? REAL(x00*y00) : REAL(alpha*x00*y00);
	    if (add) *A.ptr() += temp;
	    else *A.ptr() = temp;
	  } else {
	    T temp = a1 ? (x00 * y00) : (alpha * (x00 * y00));
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

	  Matrix<T> tempA10 = x10 * y00;
	  tempA10 += x11 * y10;

	  RecursiveInPlaceSymMultMM<ha,a1,add>(alpha,x11,y11,A11);
	  RecursiveSymMultMM<ha,a1,true>(alpha,x10,y01,A11);
	  RecursiveInPlaceSymMultMM<ha,a1,add>(alpha,x00,y00,A00);
	  RecursiveSymMultMM<ha,a1,true>(alpha,x01,y10,A00);

	  if (add) A10 += alpha * tempA10;
	  else A10 = alpha * tempA10;

	}
      }
							    
							    
  template <class T, class Tx, class Ty> inline void InPlaceSymMultMM(
      const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
      const int beta, const SymMatrixView<T>& A)
  {
    if (A.isherm())
      if (alpha == T(1))
	if (beta == 0)
	  RecursiveInPlaceSymMultMM<true,true,false>(alpha,x,y,A); 
	else
	  RecursiveInPlaceSymMultMM<true,true,true>(alpha,x,y,A); 
      else
	if (beta == 0)
	  RecursiveInPlaceSymMultMM<true,false,false>(alpha,x,y,A); 
	else
	  RecursiveInPlaceSymMultMM<true,false,true>(alpha,x,y,A); 
    else
      if (alpha == T(1))
	if (beta == 0)
	  RecursiveInPlaceSymMultMM<false,true,false>(alpha,x,y,A); 
	else
	  RecursiveInPlaceSymMultMM<false,true,true>(alpha,x,y,A); 
      else
	if (beta == 0)
	  RecursiveInPlaceSymMultMM<false,false,false>(alpha,x,y,A); 
	else
	  RecursiveInPlaceSymMultMM<false,false,true>(alpha,x,y,A); 
  }

  template <class T, class Tx, class Ty> inline void DoSymMultMM(
      const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
      const int beta, const SymMatrixView<T>& A)
  { 
    if (x.Real().cptr() == A.Real().cptr() || 
	y.Real().cptr() == A.Real().cptr()) {
      const size_t N = A.size();
      TMVAssert(x.colsize() == N);
      TMVAssert(y.rowsize() == N);
      const size_t K = x.rowsize();
      TMVAssert(y.colsize() == K);

      if (K >= N) {
	InPlaceSymMultMM(alpha,x.Cols(0,N),y.Rows(0,N),beta,A);
	if (K > N) DoSymMultMM(alpha,x.Cols(N,K),y.Rows(N,K),beta,A);
      } else { 
	DoSymMultMM(alpha,x.Rows(K,N),y.Cols(K,N),beta,
	    A.SubSymMatrix(K,N));
	if (y.stepi() < y.stepj() != A.stepi() < A.stepj()) {
	  if (x.stepi() < x.stepj() == A.stepi() < A.stepj()) {
	    // Then both x and y overlap with A(K:N,0:K)
	    // Need a temporary
	    if (A.iscm()) {
	      Matrix<T,ColMajor> temp = x.Rows(K,N) * y.Cols(0,K);
	      AddMM(alpha,temp,T(beta),A.SubMatrix(K,N,0,K));
	    } else {
	      Matrix<T,RowMajor> temp = x.Rows(K,N) * y.Cols(0,K);
	      AddMM(alpha,temp,T(beta),A.SubMatrix(K,N,0,K));
	    }
	  } else {
	    MultMM(alpha,x.SubMatrix(K,N,K,N),y.SubMatrix(K,N,0,K),T(beta),
		A.SubMatrix(K,N,0,K) );
	    MultMM(alpha,x.SubMatrix(K,N,0,K),y.SubMatrix(0,K,0,K),T(1),
		A.SubMatrix(K,N,0,K) );
	  }
	}
	else {
	  MultMM(alpha,x.Rows(K,N),y.Cols(0,K),T(beta),A.SubMatrix(K,N,0,K) );
	}
	InPlaceSymMultMM(alpha,x.Rows(0,K),y.Rows(0,K),beta,
	    A.SubSymMatrix(0,K));
      }
    }
    else 
      if (A.isherm())
	if (alpha == T(1))
	  if (beta == 0)
	    RecursiveSymMultMM<true,true,false>(alpha,x,y,A); 
	  else
	    RecursiveSymMultMM<true,true,true>(alpha,x,y,A); 
	else
	  if (beta == 0)
	    RecursiveSymMultMM<true,false,false>(alpha,x,y,A); 
	  else
	    RecursiveSymMultMM<true,false,true>(alpha,x,y,A); 
      else
	if (alpha == T(1))
	  if (beta == 0)
	    RecursiveSymMultMM<false,true,false>(alpha,x,y,A); 
	  else
	    RecursiveSymMultMM<false,true,true>(alpha,x,y,A); 
	else
	  if (beta == 0)
	    RecursiveSymMultMM<false,false,false>(alpha,x,y,A); 
	  else
	    RecursiveSymMultMM<false,false,true>(alpha,x,y,A); 
  }

  template <class T, class Tx, class Ty> void SymMultMM(
      const T alpha, const GenMatrix<Tx>& x, const GenMatrix<Ty>& y,
      const int beta, const SymMatrixView<T>& A)
    // A = beta*A + alpha * x * y
  {
    TMVAssert(A.size() == x.colsize());
    TMVAssert(A.size() == y.rowsize());
    TMVAssert(x.rowsize() == y.colsize());

#ifdef XDEBUG
    Matrix<T> A0 = A;
    Matrix<Tx> x0 = x;
    Matrix<Ty> y0 = y;
    Matrix<T> A2(A);
    if (beta == 0) A2 = alpha*x*y;
    else A2 += alpha*x*y;
#endif

    TMVAssert(beta == 0 || beta == 1);
    if (alpha != T(0) && A.size() > 0) {
      if (A.uplo() == Upper) {
	if (A.isherm()) SymMultMM(alpha,x,y,beta,A.Adjoint());
	else SymMultMM(alpha,x,y,beta,A.Transpose());
      }
      else if (A.isconj()) 
	SymMultMM(CONJ(alpha),x.Conjugate(),y.Conjugate(),beta,
	    A.Conjugate());
      else DoSymMultMM(alpha,x,y,beta,A);
    }
  
#ifdef XDEBUG
    if (Norm(A-A2) > 0.001*(abs(alpha)*Norm(x0)*Norm(y0)+
	  (beta==0?RealType(T)(0):Norm(A0)))) {
      cerr<<"SymMultMM: alpha,beta = "<<alpha<<','<<beta<<endl;
      cerr<<"x = "<<Type(x)<<"  "<<x0<<endl;
      cerr<<"y = "<<Type(y)<<"  "<<y0<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"x*y = "<<x0*y0<<endl;
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


