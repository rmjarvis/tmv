
#include "TMV_Herm.h"

namespace tmv {

  //
  // MultMV
  //
  template <class T1, class T2, class IT2, class T3, class IT3, class T4, class T5, class IT5> 
    void NonBlasMultMV(const T1 alpha, const GenHermMatrix<T2,IT2>& A, 
        const GenVector<T3,IT3>& x, 
        const T4 beta, const ModSubVector<T5,IT5>& y,
	const size_t start=0, const size_t end=A.size()) 
  // y = alpha * A * x + beta * y
  {
    TMVAssert(start<end);
    TMVAssert(end<=A.size());
    TMVAssert(end-start == y.size());
    TMVAssert(x.size() == y.size());
    typename IT5::VIT yit = y.begin();
    
    if (beta == T4(0)) 
      if (alpha == T1(1)) 
	for(size_t i=start;i<end;++i,++yit) 
	  *yit = A.row_a(i)*x.SubVector(0,i) +
	    A.row_b(i)*x.SubVector(i,x.size());
      else 
	for(size_t i=start;i<end;++i,++yit) 
	  *yit = alpha*(A.row_a(i)*x.SubVector(0,i) +
	    A.row_b(i)*x.SubVector(i,x.size()));
    else if (beta == T4(1)) 
      if (alpha == T1(1)) 
	for(size_t i=start;i<end;++i,++yit) 
	  *yit += A.row_a(i)*x.SubVector(0,i) +
	    A.row_b(i)*x.SubVector(i,x.size());
      else 
	for(size_t i=start;i<end;++i,++yit) 
	  *yit += alpha*(A.row_a(i)*x.SubVector(0,i) +
	    A.row_b(i)*x.SubVector(i,x.size()));
    else 
      if (alpha == T1(1)) 
	for(size_t i=start;i<end;++i,++yit) 
	  *yit = A.row_a(i)*x.SubVector(0,i) +
	    A.row_b(i)*x.SubVector(i,x.size()) + beta * (*yit);
      else
	for(size_t i=start;i<end;++i,++yit) 
	  *yit = alpha*(A.row_a(i)*x.SubVector(0,i) +
	    A.row_b(i)*x.SubVector(i,x.size())) + beta * (*yit);
  }

  template <class T1, class T2, class IT2, class T3, class T4, class IT4>
    void AddMM(const T1 alpha, const GenHermMatrix<T2,IT2>& A,
	const T3 beta, const ModSubMatrix<T4,IT4>& B)
    // B = alpha * A + beta * B
    {
      TMVAssert(A.colsize() == B.colsize());
      TMVAssert(A.rowsize() == B.rowsize());
      if (beta == T3(0)) {
	if (alpha == T1(1)) 
	  for(size_t i=0;i<A.colsize();++i) {
	    B.row(i).SubVector(0,i) = A.row_a(i);
	    B.row(i).SubVector(i,size()) = A.row_b(i);
	  }
	else if (alpha == T1(0)) B.Zero();
	else 
	  for(size_t i=0;i<A.colsize();++i) {
	    B.row(i).SubVector(0,i) = alpha * A.row_a(i);
	    B.row(i).SubVector(i,size()) = alpha * A.row_b(i);
	  }
      } else {
	if (beta != T3(1)) B *= beta;
	if (alpha == T1(1))
	  for(size_t i=0;i<A.colsize();++i) {
	    B.row(i).SubVector(0,i) += A.row_a(i);
	    B.row(i).SubVector(i,size()) += A.row_b(i);
	  }
	else
	  for(size_t i=0;i<A.colsize();++i) {
	    B.row(i).SubVector(0,i) += alpha*A.row_a(i);
	    B.row(i).SubVector(i,size()) += alpha*A.row_b(i);
	  }
      }
    }

  template <class T1, class T2, class IT2, class T3, class T4, class IT4>
    void AddMM(const T1 alpha, const GenHermMatrix<T2,IT2>& A,
	const T3 beta, const ModSubHermMatrix<T4,IT4>& B)
    // B = alpha * A + beta * B
    {
      TMVAssert(A.colsize() == B.colsize());
      TMVAssert(A.rowsize() == B.rowsize());
      if (beta == T3(0)) {
	if (alpha == T1(1)) 
	  for(size_t i=0;i<A.colsize();++i) B.row_a(i) = A.row_a(i);
	else if (alpha == T1(0)) B.Zero();
	else 
	  for(size_t i=0;i<A.colsize();++i) B.row_a(i) = alpha*A.row_a(i);
      } else if (beta == T3(1)) {
	if (alpha == T1(1))
	  for(size_t i=0;i<A.colsize();++i) B.row_a(i) += A.row_a(i);
	else if (alpha != T1(0)) 
	  for(size_t i=0;i<A.colsize();++i) B.row_a(i) += alpha*A.row_a(i);
      } else {
	if (alpha == T1(1))
	  for(size_t i=0;i<A.colsize();++i) B.row_a(i) = A.row_a(i) +
	    beta*B.row_a(i);;
	else if (alpha == T1(0))
	  for(size_t i=0;i<A.colsize();++i) B.row_a(i) *= beta;
	else
	  for(size_t i=0;i<A.colsize();++i) B.row_a(i) = alpha*A.row_a(i) +
	    beta*B.row_a(i);;
      }
    }

  template <class T1, class T2, class IT2, class T3, class IT3, class T4, class T5, class IT5> 
    void NonBlasMultMM(const T1 alpha,
	const GenHermMatrix<T2,IT2>& A, const GenMatrix<T3,IT3>& B, 
	const T4 beta, const ModSubMatrix<T5,IT5>& C) 
    {
      TMVAssert(A.rowsize() == B.colsize());
      TMVAssert(A.colsize() == C.colsize());
      TMVAssert(B.rowsize() == C.rowsize());
      for(size_t i=0;i<C.rowsize();++i)
	MultMV(alpha,A,B.col(i),beta,C.col(i));
    }
  template <class T1, class T2, class IT2, class T3, class IT3, class T4, class T5, class IT5> 
    void NonBlasMultMM(const T1 alpha,
        const GenMatrix<T2,IT2>& A, const GenHermMatrix<T3,IT3>& B, 
	const T4 beta, const ModSubMatrix<T5,IT5>& C) 
    {
      TMVAssert(A.rowsize() == B.colsize());
      TMVAssert(A.colsize() == C.colsize());
      TMVAssert(B.rowsize() == C.rowsize());
      for(size_t i=0;i<C.rowsize();++i)
	MultMV(alpha,B.Transpose(),A.row(i),beta,C.row(i));
    }
#ifdef BLAS
  template <class T1, class T2, class IT2, class T3, class IT3, class T4, class T5, class IT5> 
    void BlasMultMM(const T1 alpha, 
	const GenHermMatrix<T2,IT2>& A, const GenMatrix<T3,IT3>& B,
	const T4 beta, const ModSubMatrix<T5,IT5>& C) 
    { NonBlasMultMM(alpha,A,B,beta,C); }
  template <class T1, class T2, class IT2, class T3, class IT3, class T4, class T5, class IT5> 
    void BlasMultMM(const T1 alpha, 
	const GenMatrix<T2,IT2>& A, const GenHermMatrix<T3,IT3>& B,
	const T4 beta, const ModSubMatrix<T5,IT5>& C) 
    { NonBlasMultMM(alpha,A,B,beta,C); }
  template <> void BlasMultMM(const double alpha,
      const GenHermMatrix<double>& A, const GenMatrix<double>& B,
      const double beta, const ModSubMatrix<double>& C) 
  { 
    cblas_dsymm(CblasRowMajor,CblasLeft,CblasLower,
	C.colsize(),C.rowsize(),alpha,A.cptr(),A.stepi(),
	B.cptr(),B.stepi(),beta,C.ptr(),C.stepi());
  }
  template <> void BlasMultMM(const double alpha,
      const GenMatrix<double>& A, const GenHermMatrix<double>& B,
      const double beta, const ModSubMatrix<double>& C) 
  { 
    cblas_dsymm(CblasRowMajor,CblasRight,CblasLower,
	C.colsize(),C.rowsize(),alpha,B.cptr(),B.stepi(),
	A.cptr(),A.stepi(),beta,C.ptr(),C.stepi());
  }
  template <class T1, class T4, class IT> void BlasMultMM(const T1 alpha, 
      const GenHermMatrix<complex<double>,IT>& A,
      const GenMatrix<complex<double>,IT>& B, 
      const T4 beta, const ModSubMatrix<complex<double>,IT>& C) 
  { 
    complex<double> calpha = alpha;
    complex<double> cbeta = beta;
    cblas_zhemm(CblasRowMajor,CBlasLeft,CblasLower,
	C.colsize(),C.rowsize(),&calpha,A.cptr(),A.stepi(),
	B.cptr(),B.stepi(),&cbeta,C.ptr(),C.stepi());
  }
  template <class T1, class T4, class IT> void BlasMultMM(const T1 alpha, 
      const GenMatrix<complex<double>,IT>& A, 
      const GenHermMatrix<complex<double>,IT>& B,
      const T4 beta, const ModSubMatrix<complex<double>,IT>& C) 
  { 
    complex<double> calpha = alpha;
    complex<double> cbeta = beta;
    cblas_zhemm(CblasRowMajor,CBlasLeft,CblasLower,
	C.colsize(),C.rowsize(),&calpha,B.cptr(),B.stepi(),
	A.cptr(),A.stepi(),&cbeta,C.ptr(),C.stepi());
  }
#ifndef NOFLOAT
  template <> void BlasMultMM(const float alpha,
      const GenHermMatrix<float>& A, const GenMatrix<float>& B,
      const float beta, const ModSubMatrix<float>& C) 
  { 
    cblas_dsymm(CblasRowMajor,CblasLeft,CblasLower,
	C.colsize(),C.rowsize(),alpha,A.cptr(),A.stepi(),
	B.cptr(),B.stepi(),beta,C.ptr(),C.stepi());
  }
  template <> void BlasMultMM(const float alpha,
      const GenMatrix<float>& A, const GenHermMatrix<float>& B,
      const float beta, const ModSubMatrix<float>& C) 
  { 
    cblas_ssymm(CblasRowMajor,CblasRight,CblasLower,
	C.colsize(),C.rowsize(),alpha,B.cptr(),B.stepi(),
	A.cptr(),A.stepi(),beta,C.ptr(),C.stepi());
  }
  template <class T1, class T4, class IT> void BlasMultMM(const T1 alpha, 
      const GenHermMatrix<complex<float>,IT>& A,
      const GenMatrix<complex<float>,IT>& B, 
      const T4 beta, const ModSubMatrix<complex<float>,IT>& C) 
  { 
    complex<float> calpha = alpha;
    complex<float> cbeta = beta;
    cblas_chemm(CblasRowMajor,CBlasLeft,CblasLower,
	C.colsize(),C.rowsize(),&calpha,A.cptr(),A.stepi(),
	B.cptr(),B.stepi(),&cbeta,C.ptr(),C.stepi());
  }
  template <class T1, class T4, class IT> void BlasMultMM(const T1 alpha, 
      const GenMatrix<complex<float>,IT>& A, 
      const GenHermMatrix<complex<float>,IT>& B,
      const T4 beta, const ModSubMatrix<complex<float>,IT>& C) 
  { 
    complex<float> calpha = alpha;
    complex<float> cbeta = beta;
    cblas_chemm(CblasRowMajor,CBlasLeft,CblasLower,
	C.colsize(),C.rowsize(),&calpha,B.cptr(),B.stepi(),
	A.cptr(),A.stepi(),&cbeta,C.ptr(),C.stepi());
  }
#endif
#endif // BLAS

  template <class T1, class T2, class IT2, class T3, class IT3, class T4, class T5, class IT5> 
    void MultMM(const T1 alpha, const GenHermMatrix<T2,IT2>& A, 
        const GenMatrix<T3,IT3>& B, 
	const T4 beta, const ModSubMatrix<T5,IT5>& C) 
    // y = alpha * A * B + beta * C
    {
      TMVAssert(A.rowsize() == B.colsize());
      TMVAssert(A.colsize() == C.colsize());
      TMVAssert(B.rowsize() == C.rowsize());
#ifdef BLAS
      if (A.stepj() == 1) {
	if (B.SameStorageAs(C)) {
	  if (B.SameAs(C)) {
	    NonBlasMultMM(alpha,A,B,beta,C);
	  } else {
	    if (C.stepi() == 1) {
	      Matrix<T3> tempB(B.Transpose());
	      BlasMultMM(alpha,tempB,A.Transpose(),beta,C.Transpose());
	    } else if (C.stepj() == 1) {
	      Matrix<T3> tempB(B);
	      BlasMultMM(alpha,A,tempB,beta,C);
	    } else {
	      NonBlasMultMM(alpha,A,B,beta,C);
	    }
	  } 
	} else {
	  if (C.stepi() == 1 && B.stepi() == 1) {
	    BlasMultMM(alpha,B.Transpose(),A.Transpose(),beta,C.Transpose());
	  } else if (C.stepj() == 1 && B.stepj() == 1) {
	    BlasMultMM(alpha,A,B,beta,C);
	  } else {
	    NonBlasMultMM(alpha,A,B,beta,C);
	  }
	}
      } else 
#endif
	NonBlasMultMM(alpha,A,B,beta,C);
    }

  template <class T1, class T2, class IT2, class T3, class IT3, class T4, class T5, class IT5> 
    void MultMM(const T1 alpha,
        const GenMatrix<T3,IT3>& A, const GenHermMatrix<T2,IT2>& B, 
	const T4 beta, const ModSubMatrix<T5,IT5>& C) 
    // y = alpha * A * B + beta * C
    {
      TMVAssert(A.rowsize() == B.colsize());
      TMVAssert(A.colsize() == C.colsize());
      TMVAssert(B.rowsize() == C.rowsize());
#ifdef BLAS
      if (B.stepj() == 1) {
	if (A.SameStorageAs(C)) {
	  if (A.SameAs(C)) {
	    NonBlasMultMM(alpha,A,B,beta,C);
	  } else {
	    if (C.stepi() == 1) {
	      Matrix<T3> tempA(A.Transpose());
	      BlasMultMM(alpha,B.Transpose(),tempA,beta,C.Transpose());
	    } else if (C.stepj() == 1) {
	      Matrix<T3> tempA(A);
	      BlasMultMM(alpha,tempA,B,beta,C);
	    } else {
	      NonBlasMultMM(alpha,A,B,beta,C);
	    }
	  } 
	} else {
	  if (C.stepi() == 1 && A.stepi() == 1) {
	    BlasMultMM(alpha,B.Transpose(),A.Transpose(),beta,C.Transpose());
	  } else if (C.stepj() == 1 && B.stepj() == 1) {
	    BlasMultMM(alpha,A,B,beta,C);
	  } else {
	    NonBlasMultMM(alpha,A,B,beta,C);
	  }
	}
      } else 
#endif
	NonBlasMultMM(alpha,A,B,beta,C);
    }

  template <class T1, class T2, class IT2, class T3, class IT3, class T4, class T5, class IT5> 
    void MultMM(const T1 alpha, const GenHermMatrix<T2,IT2>& A, 
        const GenHermMatrix<T3,IT3>& B, 
	const T4 beta, const ModSubHermMatrix<T5,IT5>& C) 
    // C = alpha * A * B + beta * C
    {
      TMVAssert(A.rowsize() == B.colsize());
      TMVAssert(A.colsize() == C.colsize());
      TMVAssert(B.rowsize() == C.rowsize());
      if (A.SameStorageAs(C)) {
	for(size_t j=0;j<C.colsize();++j) {
	  // Take the adjoint of the equation: (Use At = A, etc.)
	  // C = (alpha*) * B * A + (beta*) * C
	  MultMV(CONJ(alpha),B,A.col_b(j),CONJ(beta),C.col_b(j),j,B.colsize());
	}
      } else {
	for(size_t j=0;j<C.rowsize();++j) {
	  MultMV(alpha,A,B.col_b(j),beta,C.col_b(j),j,A.colsize());
	}
      }
    }

#define InstFile "TMV_HermMatrixArith.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


