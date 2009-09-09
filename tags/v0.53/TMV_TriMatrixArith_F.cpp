
#include "TMV.h"
#include "TMV_Tri.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
  const size_t TRI_MM_BLOCKSIZE = TMV_BLOCKSIZE;
  const size_t TRI_MM_BLOCKSIZE2 = TMV_BLOCKSIZE/2;
#else
  const size_t TRI_MM_BLOCKSIZE = 64;
  const size_t TRI_MM_BLOCKSIZE2 = 32;
#endif

  //
  // MultMM: U = U * U
  //

  template <class T, class Ta> void RRMultEqMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, const UpperTriMatrixView<T>& B)
  {
#ifdef XDEBUG
    //cerr<<"RRMultEqMM Upper/Upper alpha = "<<alpha<<endl;
    Matrix<T> B2 = alpha*Matrix<T>(A)*Matrix<T>(B);
    Matrix<T> B0 = B;
#endif
    TMVAssert(A.isrm());
    TMVAssert(B.isrm());
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
    TMVAssert(B.size() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(B.ct()==NonConj);

    const size_t N = B.size();

    if (A.isunit()) {
      for(size_t i=0; i<N; ++i) {
	B.row(i,i+1,N) += A.row(i,i+1,N) * B.SubTriMatrix(i+1,N);
	if (!B.isunit()) 
	  B.row(i,i,N) *= alpha;
	else { TMVAssert(alpha == T(1)); }
      }
    }
    else {
      TMVAssert(!B.isunit());
      const int Ads = A.stepi()+A.stepj();
      const Ta* Aii = A.cptr();
      const int Bds = B.stepi()+B.stepj();
      T* Bii = B.ptr();
      for(size_t i=0; i<N; ++i,Aii+=Ads,Bii+=Bds) {
	T aa = A.isconj()?CONJ(*Aii):*Aii;
	if (alpha != T(1)) aa *= alpha;
	B.row(i,i+1,N) = aa * B.row(i,i+1,N) +
	    alpha * A.row(i,i+1,N) * B.SubTriMatrix(i+1,N);
	*Bii *= aa;
      }
    }
#ifdef XDEBUG
    if (Norm(B-B2) > 0.001*(Norm(A)+Norm(B))) {
      cerr<<"RRMultEqMM alpha = "<<alpha<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"--> B = "<<B<<endl;
      cerr<<"B2 = "<<B2<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta> void CRMultEqMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, const UpperTriMatrixView<T>& B)
  {
#ifdef XDEBUG
    //cerr<<"CRMultEqMM Upper/Upper alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.iscm());
    TMVAssert(B.isrm());
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
    TMVAssert(B.size() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(B.ct()==NonConj);

    const size_t N = B.size();

    if (A.isunit()) {
      if (B.isunit()) {
	TMVAssert(alpha == T(1));
	for(size_t j=1; j<N; ++j) {
	  B.SubMatrix(0,j,j+1,N) += A.col(j,0,j) ^ B.row(j,j+1,N);
	  B.col(j,0,j) += A.col(j,0,j);
	}
      }
      else {
	for(size_t j=0; j<N; ++j) 
	  B.SubMatrix(0,j,j,N) += A.col(j,0,j) ^ B.row(j,j,N);
	B *= alpha;
      }
    } else {
      TMVAssert(!B.isunit());
      const int Ads = A.stepi()+A.stepj();
      const Ta* Ajj = A.cptr();
      for(size_t j=0; j<N; ++j,Ajj+=Ads) {
	B.SubMatrix(0,j,j,N) += A.col(j,0,j) ^ B.row(j,j,N);
	B.row(j,j,N) *= A.isconj()?CONJ(*Ajj):*Ajj;
      }
      B *= alpha;
    } 
  }

  template <class T, class Ta> void CMultEqMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, const UpperTriMatrixView<T>& B)
  {
#ifdef XDEBUG
    //cerr<<"CMultEqMM Upper/Upper alpha = "<<alpha<<endl;
#endif
    TMVAssert(B.iscm());
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
    TMVAssert(B.size() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(B.ct()==NonConj);

    if (B.isunit()) {
      // Then alpha = 1 and A.isunit
      for(size_t j=1;j<B.size();++j) {
	B.col(j,0,j) = A.SubTriMatrix(0,j) * B.col(j,0,j);
	B.col(j,0,j) += A.col(j,0,j);
      }
    }
    else {
      for(size_t j=0;j<B.size();++j) 
	B.col(j,0,j+1) = alpha * A.SubTriMatrix(0,j+1) * B.col(j,0,j+1);
    }
  }

  template <class T, class Ta> inline void MultEqMM(T alpha,
      const GenUpperTriMatrix<Ta>& A, const UpperTriMatrixView<T>& B)
    // B = alpha * A * B
  {
#ifdef XDEBUG
    //cerr<<"MultEqMM: "<<alpha<<"  "<<Type(A)<<"  "<<A<<"  "<<Type(B)<<"  "<<B<<endl;
#endif
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
    TMVAssert(alpha != T(0));
    TMVAssert(B.ct()==NonConj);
    TMVAssert(B.size() > 0);

    const size_t nb = TRI_MM_BLOCKSIZE;
    const size_t N = A.size();

    if (N <= TRI_MM_BLOCKSIZE2) {
      if (A.isrm() && B.isrm())
	RRMultEqMM(alpha,A,B);
      else if (A.iscm() && B.isrm())
	CRMultEqMM(alpha,A,B);
      else if (B.iscm())
	CMultEqMM(alpha,A,B);
      else {
	if (B.isunit()) {
	  UpperTriMatrix<T,UnitDiag,ColMajor> BB = B;
	  if (!(A.isrm() || A.iscm())) {
	    if (A.isunit()) {
	      UpperTriMatrix<T,UnitDiag,ColMajor> AA = A;
	      CMultEqMM(alpha,AA,BB.View());
	    } else {
	      UpperTriMatrix<T,NonUnitDiag,ColMajor> AA = A;
	      CMultEqMM(alpha,AA,BB.View());
	    }
	  }
	  else CMultEqMM(alpha,A,BB.View());
	  B = BB;
	} else {
	  UpperTriMatrix<T,NonUnitDiag,ColMajor> BB = B;
	  if (!(A.isrm() || A.iscm())) {
	    if (A.isunit()) {
	      UpperTriMatrix<T,UnitDiag,ColMajor> AA = A;
	      CMultEqMM(alpha,AA,BB.View());
	    } else {
	      UpperTriMatrix<T,NonUnitDiag,ColMajor> AA = A;
	      CMultEqMM(alpha,AA,BB.View());
	    }
	  }
	  else CMultEqMM(alpha,A,BB.View());
	  B = BB;
	}
      }
    } else {
      size_t k = N/2;
      if (k > nb) k = k/nb*nb;

      // [ A00 A01 ] [ B00 B01 ] = [ A00 B00   A00 B01 + A01 B11 ]
      // [  0  A11 ] [  0  B11 ]   [    0           A11 B11      ]

      ConstUpperTriMatrixView<Ta> A00 = A.SubTriMatrix(0,k);
      ConstMatrixView<Ta> A01 = A.SubMatrix(0,k,k,N);
      ConstUpperTriMatrixView<Ta> A11 = A.SubTriMatrix(k,N);
      UpperTriMatrixView<T> B00 = B.SubTriMatrix(0,k);
      MatrixView<T> B01 = B.SubMatrix(0,k,k,N);
      UpperTriMatrixView<T> B11 = B.SubTriMatrix(k,N);

      MultEqMM(alpha,A00,B00);
      B01 = alpha * A00 * B01;
      B01 += alpha * A01 * B11;
      MultEqMM(alpha,A11,B11);
    }
  }

  template <class T, class Ta> void RRMultEqMM(
      T alpha, const GenLowerTriMatrix<Ta>& A, const LowerTriMatrixView<T>& B)
  {
#ifdef XDEBUG
    //cerr<<"RRMultEqMM Lower/Lower: alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.isrm());
    TMVAssert(B.isrm());
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
    TMVAssert(B.size() > 0);
    TMVAssert(alpha != T(0));

    const size_t N = B.size();

    if (A.isunit()) {
      for(int i=N-1; i>=0; --i) {
	B.row(i,0,i) += A.row(i,0,i) * B.SubTriMatrix(0,i);
	B.row(i,0,i) *= alpha;
      }
      if (!B.isunit()) 
	B.diag() *= alpha;
    }
    else {
      TMVAssert(!B.isunit());
      const int Ads = A.stepi()+A.stepj();
      const Ta* Aii = A.cptr()+(N-1)*Ads;
      const int Bds = B.stepi()+B.stepj();
      T* Bii = B.ptr()+(N-1)*Bds;
      for(size_t i=N-1; i>0; --i,Aii-=Ads,Bii-=Bds) {
	T aa = A.isconj()?CONJ(*Aii):*Aii;
	if (alpha != T(1)) aa *= alpha;
	B.row(i,0,i) = aa * B.row(i,0,i) +
	  alpha * A.row(i,0,i) * B.SubTriMatrix(0,i);
	*Bii *= aa;
      }
    }
  }

  template <class T, class Ta> void CRMultEqMM(
      T alpha, const GenLowerTriMatrix<Ta>& A, const LowerTriMatrixView<T>& B)
  {
#ifdef XDEBUG
    //cerr<<"CRMultEqMM Lower/Lower: alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.iscm());
    TMVAssert(B.isrm());
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
    TMVAssert(B.size() > 0);
    TMVAssert(alpha != T(0));

    const size_t N = B.size();

    if (A.isunit()) {
      if (B.isunit()) {
	for(int j=N-1; j>=0; --j) {
	  B.SubMatrix(j+1,N,0,j) += A.col(j,j+1,N) ^ B.row(j,0,j);
	  B.col(j,j+1,N) += A.col(j,j+1,N);
	}
      }
      else {
	for(int j=N-1; j>=0; --j) 
	  B.SubMatrix(j+1,N,0,j+1) += A.col(j,j+1,N) ^ B.row(j,0,j+1);
      }
    } else {
      const int Ads = A.stepi()+A.stepj();
      const Ta* Ajj = A.cptr()+(N-1)*Ads;
      TMVAssert(!B.isunit());
      for(int j=N-1; j>=0; --j,Ajj-=Ads) {
	B.SubMatrix(j+1,N,0,j+1) += A.col(j,j+1,N) ^ B.row(j,0,j+1);
	B.row(j,0,j+1) *= A.isconj()?CONJ(*Ajj):*Ajj;
      }
    } 
    B *= alpha;
  }

  template <class T, class Ta> void CMultEqMM(
      T alpha, const GenLowerTriMatrix<Ta>& A, const LowerTriMatrixView<T>& B)
  {
#ifdef XDEBUG
    //cerr<<"CMultEqMM Lower/Lower: alpha = "<<alpha<<endl;
#endif
    TMVAssert(B.iscm());
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
    TMVAssert(B.size() > 0);
    TMVAssert(alpha != T(0));

    const size_t N = B.size();
    if (B.isunit()) {
      // Then alpha = 1 and A.isunit
      for(size_t j=0;j<B.size()-1;++j) {
	B.col(j,j+1,N) = A.SubTriMatrix(j+1,N) * B.col(j,j+1,N);
	B.col(j,j+1,N) += A.col(j,j+1,N);
      }
    }
    else {
      for(size_t j=0;j<B.size();++j) 
	B.col(j,j,N) = alpha * A.SubTriMatrix(j,N) * B.col(j,j,N);
    }
  }

  template <class T, class Ta> inline void MultEqMM(T alpha,
      const GenLowerTriMatrix<Ta>& A, const LowerTriMatrixView<T>& B)
    // B = alpha * A * B
  {
#ifdef XDEBUG
    //cerr<<"MultEqMM: "<<alpha<<"  "<<Type(A)<<"  "<<A<<"  "<<Type(B)<<"  "<<B<<endl;
#endif
    TMVAssert(A.size() == B.size());
    TMVAssert(!B.isunit() || (A.isunit() && alpha == T(1)) );
    TMVAssert(alpha != T(0));
    TMVAssert(B.ct()==NonConj);
    TMVAssert(B.size() > 0);

    const size_t nb = TRI_MM_BLOCKSIZE;
    const size_t N = A.size();

    if (N <= TRI_MM_BLOCKSIZE2) {
      if (A.isrm() && B.isrm())
	RRMultEqMM(alpha,A,B);
      else if (A.iscm() && B.isrm())
	CRMultEqMM(alpha,A,B);
      else if (B.iscm())
	CMultEqMM(alpha,A,B);
      else {
	if (B.isunit()) {
	  LowerTriMatrix<T,UnitDiag,ColMajor> BB = B;
	  if (!(A.isrm() || A.iscm())) {
	    if (A.isunit()) {
	      LowerTriMatrix<T,UnitDiag,ColMajor> AA = A;
	      CMultEqMM(alpha,AA,BB.View());
	    } else {
	      LowerTriMatrix<T,NonUnitDiag,ColMajor> AA = A;
	      CMultEqMM(alpha,AA,BB.View());
	    }
	  }
	  else CMultEqMM(alpha,A,BB.View());
	  B = BB;
	} else {
	  LowerTriMatrix<T,NonUnitDiag,ColMajor> BB = B;
	  if (!(A.isrm() || A.iscm())) {
	    if (A.isunit()) {
	      LowerTriMatrix<T,UnitDiag,ColMajor> AA = A;
	      CMultEqMM(alpha,AA,BB.View());
	    } else {
	      LowerTriMatrix<T,NonUnitDiag,ColMajor> AA = A;
	      CMultEqMM(alpha,AA,BB.View());
	    }
	  }
	  else CMultEqMM(alpha,A,BB.View());
	  B = BB;
	}
      }
    } else {
      size_t k = N/2;
      if (k > nb) k = k/nb*nb;

      // [ A00  0  ] [ B00  0  ] = [      A00 B00           0    ]
      // [ A10 A11 ] [ B10 B11 ]   [ A10 B00 + A11 B10   A11 B11 ]

      ConstLowerTriMatrixView<Ta> A00 = A.SubTriMatrix(0,k);
      ConstMatrixView<Ta> A10 = A.SubMatrix(k,N,0,k);
      ConstLowerTriMatrixView<Ta> A11 = A.SubTriMatrix(k,N);
      LowerTriMatrixView<T> B00 = B.SubTriMatrix(0,k);
      MatrixView<T> B10 = B.SubMatrix(k,N,0,k);
      LowerTriMatrixView<T> B11 = B.SubTriMatrix(k,N);

      MultEqMM(alpha,A11,B11);
      B10 = alpha * A11 * B10;
      B10 += alpha * A10 * B00;
      MultEqMM(alpha,A00,B00);
    }
  }

  template <class T, class Ta, class Tb> void RRAddMultMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, 
      const GenUpperTriMatrix<Tb>& B, const UpperTriMatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"RRAddMultMM Upper/Upper\n";
    Matrix<T> C2 = Matrix<T>(C)+alpha*Matrix<T>(A)*Matrix<T>(B);
    Matrix<T> C0 = C;
#endif
    TMVAssert(A.isrm());
    TMVAssert(B.isrm());
    TMVAssert(C.isrm());
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.size());
    TMVAssert(!C.isunit());
    TMVAssert(C.size() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct()==NonConj);

    const size_t N = C.size();

    if (A.isunit()) {
      if (B.isunit()) {
	const int Cds = C.stepi()+C.stepj();
	T* Cii = C.ptr();
	for(size_t i=0; i<N; ++i,Cii+=Cds) {
	  C.row(i,i+1,N) += alpha * A.row(i,i+1,N) * B.SubTriMatrix(i+1,N);
	  C.row(i,i+1,N) += alpha * B.row(i,i+1,N);
	  *Cii += alpha;
	}
      } else {
	for(size_t i=0; i<N; ++i) {
	  C.row(i,i+1,N) += alpha * A.row(i,i+1,N) * B.SubTriMatrix(i+1,N);
	  C.row(i,i,N) += alpha * B.row(i,i,N);
	}
      }
    } else {
      for(size_t i=0; i<N; ++i) 
	C.row(i,i,N) += alpha * A.row(i,i,N) * B.SubTriMatrix(i,N);
    }
#ifdef XDEBUG
    if (Norm(C-C2) > 0.001*(Norm(C)+Norm(A)+Norm(B))) {
      cerr<<"RRAddMultMM alpha = "<<alpha<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"--> C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta, class Tb> void CRAddMultMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, const GenUpperTriMatrix<Tb>& B,
      const UpperTriMatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"CRAddMult Upper/Upper: alpha = "<<alpha<<endl;
#endif
    TMVAssert(A.iscm());
    TMVAssert(B.isrm());
    TMVAssert(C.isrm());
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.size());
    TMVAssert(!C.isunit());
    TMVAssert(C.size() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct()==NonConj);

    const size_t N = C.size();

    if (A.isunit()) {
      if (B.isunit()) {
	const int Cds = C.stepi()+C.stepj();
	T* Cjj = C.ptr();
	for(size_t j=0; j<N; ++j,Cjj+=Cds) {
	  C.SubMatrix(0,j,j+1,N) += alpha * A.col(j,0,j) ^ B.row(j,j+1,N);
	  C.col(j,0,j) += alpha * A.col(j,0,j);
	  C.row(j,j+1,N) += alpha * B.row(j,j+1,N);
	  *Cjj += alpha;
	}
      }
      else {
	for(size_t j=0; j<N; ++j) {
	  C.SubMatrix(0,j,j,N) += alpha * A.col(j,0,j) ^ B.row(j,j,N);
	  C.row(j,j,N) += alpha * B.row(j,j,N);
	}
      }
    } else {
      if (B.isunit()) {
	for(size_t j=0; j<N; ++j) {
	  C.SubMatrix(0,j+1,j+1,N) += alpha * A.col(j,0,j+1)^B.row(j,j+1,N);
	  C.row(j,0,j+1) += alpha * A.col(j,0,j+1);
	}
      } else {
	for(size_t j=0; j<N; ++j) 
	  C.SubMatrix(0,j+1,j,N) += alpha * A.col(j,0,j+1) ^ B.row(j,j,N);
      }
    } 
  }

  template <class T, class Ta, class Tb> void CAddMultMM(
      T alpha, const GenUpperTriMatrix<Ta>& A, const GenUpperTriMatrix<Tb>& B,
      const UpperTriMatrixView<T>& C)
  {
#ifdef XDEBUG
    //cerr<<"CAddMult Upper/Upper: alpha = "<<alpha<<endl;
#endif
    TMVAssert(B.iscm());
    TMVAssert(C.iscm());
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.size());
    TMVAssert(!C.isunit());
    TMVAssert(C.size() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct()==NonConj);

    const size_t N = C.size();

    if (B.isunit()) {
      if (A.isunit()) {
	const int Cds = C.stepi()+C.stepj();
	T* Cjj = C.ptr();
	for(size_t j=0;j<N;++j,Cjj+=Cds) {
	  C.col(j,0,j) += alpha * A.SubTriMatrix(0,j) * B.col(j,0,j);
	  C.col(j,0,j) += alpha * A.col(j,0,j);
	  *Cjj += alpha;
	}
      } else {
	for(size_t j=0;j<N;++j) {
	  C.col(j,0,j) += alpha * A.SubTriMatrix(0,j) * B.col(j,0,j);
	  C.col(j,0,j+1) += alpha * A.col(j,0,j+1);
	}
      }
    } else {
      for(size_t j=0;j<N;++j)
	C.col(j,0,j+1) += alpha * A.SubTriMatrix(0,j+1) * B.col(j,0,j+1);
    }
  }

  template <class T, class Ta, class Tb> inline void AddMultMM(
      T alpha, const GenUpperTriMatrix<Ta>& A,
      const GenUpperTriMatrix<Tb>& B, const UpperTriMatrixView<T>& C)
    // C += alpha * A * B
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.size());
    TMVAssert(!C.isunit());
    TMVAssert(C.size() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct()==NonConj);

    const size_t nb = TRI_MM_BLOCKSIZE;
    const size_t N = A.size();

    if (N <= TRI_MM_BLOCKSIZE2) {
      if (A.isrm() && B.isrm() && C.isrm())
	RRAddMultMM(alpha,A,B,C);
      else if (A.iscm() && B.isrm() && C.isrm())
	CRAddMultMM(alpha,A,B,C);
      else if (!A.isrm() && !A.iscm()) {
	if (A.isunit()) {
	  UpperTriMatrix<T,UnitDiag,ColMajor> AA = A;
	  AddMultMM(alpha,AA,B,C);
	} else {
	  UpperTriMatrix<T,NonUnitDiag,ColMajor> AA = A;
	  AddMultMM(alpha,AA,B,C);
	}
      }
      else if (B.iscm() && C.iscm())
	CAddMultMM(alpha,A,B,C);
      else if (C.isrm()) { 
	TMVAssert(!B.isrm());
	if (B.isunit()) {
	  UpperTriMatrix<T,UnitDiag,RowMajor> BB = B;
	  AddMultMM(alpha,A,BB,C);
	} else {
	  UpperTriMatrix<T,NonUnitDiag,RowMajor> BB = B;
	  AddMultMM(alpha,A,BB,C);
	}
      } 
      else if (C.iscm()) { 
	TMVAssert(!B.iscm());
	if (B.isunit()) {
	  UpperTriMatrix<T,UnitDiag,ColMajor> BB = B;
	  AddMultMM(alpha,A,BB,C);
	} else {
	  UpperTriMatrix<T,NonUnitDiag,ColMajor> BB = B;
	  AddMultMM(alpha,A,BB,C);
	}
      }
      else {
	UpperTriMatrix<T,NonUnitDiag,ColMajor> CC = C;
	AddMultMM(alpha,A,B,CC.View());
	C = CC;
      }
    } else {
      size_t k = N/2;
      if (k > nb) k = k/nb*nb;

      // [ A00 A01 ] [ B00 B01 ] = [ A00 B00   A00 B01 + A01 B11 ]
      // [  0  A11 ] [  0  B11 ]   [    0           A11 B11      ]

      ConstUpperTriMatrixView<Ta> A00 = A.SubTriMatrix(0,k);
      ConstMatrixView<Ta> A01 = A.SubMatrix(0,k,k,N);
      ConstUpperTriMatrixView<Ta> A11 = A.SubTriMatrix(k,N);
      ConstUpperTriMatrixView<Tb> B00 = B.SubTriMatrix(0,k);
      ConstMatrixView<Tb> B01 = B.SubMatrix(0,k,k,N);
      ConstUpperTriMatrixView<Tb> B11 = B.SubTriMatrix(k,N);
      UpperTriMatrixView<T> C00 = C.SubTriMatrix(0,k);
      MatrixView<T> C01 = C.SubMatrix(0,k,k,N);
      UpperTriMatrixView<T> C11 = C.SubTriMatrix(k,N);

      AddMultMM(alpha,A00,B00,C00);
      C01 += alpha * A00 * B01;
      C01 += alpha * A01 * B11;
      AddMultMM(alpha,A11,B11,C11);
    }
  }
  
  template <class T, class Ta, class Tb> void TempMultMM(const T alpha,
      const GenUpperTriMatrix<Ta>& A, const GenUpperTriMatrix<Tb>& B,
      const T beta, const UpperTriMatrixView<T>& C)
    // C = alpha * A * B + beta * C    
  { 
    if (B.isrm()) {
      if (B.dt() == UnitDiag) {
	UpperTriMatrix<T,UnitDiag,RowMajor> tempB(B);
	MultEqMM(alpha,A,tempB.View());
	C *= beta;
	C += tempB;
      } else {
	UpperTriMatrix<T,NonUnitDiag,RowMajor> tempB(B);
	MultEqMM(alpha,A,tempB.View());
	C *= beta;
	C += tempB;
      }
    } else {
      if (B.dt() == UnitDiag) {
	UpperTriMatrix<T,UnitDiag,ColMajor> tempB(B);
	MultEqMM(alpha,A,tempB.View());
	C *= beta;
	C += tempB;
      } else {
	UpperTriMatrix<T,NonUnitDiag,ColMajor> tempB(B);
	MultEqMM(alpha,A,tempB.View());
	C *= beta;
	C += tempB;
      }
    }
  }

  // MJ: I think the storage requirements can be relaxed here
  template <class T, class Ta, class Tb> void MultMM(const T alpha,
      const GenUpperTriMatrix<Ta>& A, const GenUpperTriMatrix<Tb>& B,
      const T beta, const UpperTriMatrixView<T>& C)
    // C = alpha * A * B + beta * C    
  { 
#ifdef XDEBUG
    //cerr<<"MultMM: "<<alpha<<"  "<<Type(A)<<"  "<<A<<"  "<<Type(B)<<"  "<<B<<endl;
    Matrix<T> C2 = beta*Matrix<T>(C) + alpha*Matrix<T>(A)*Matrix<T>(B);
    Matrix<T> C0 = C;
#endif
    TMVAssert(A.size() == C.size());
    TMVAssert(A.size() == B.size());
    TMVAssert(!C.isunit() || (A.isunit() && B.isunit() && alpha == T(1)) );

    if (C.size() > 0) {
      if (C.isconj()) MultMM(CONJ(alpha),A.Conjugate(),B.Conjugate(),
	  CONJ(beta),C.Conjugate());
      else if (alpha==T(0)) 
	C *= beta;
      else if (C.SameStorageAs(A)) 
	if (beta == T(0) && !C.SameStorageAs(B)) {
	  C = A;
	  MultEqMM(alpha,B.Transpose(),C.Transpose());
	}
	else 
	  TempMultMM(alpha,A,B,beta,C);
      else if (beta == T(0)) {
	C = B;
	MultEqMM(alpha,A,C);
      }
      else if (C.SameStorageAs(B)) 
	TempMultMM(alpha,A,B,beta,C);
      else {
	C *= beta;
	AddMultMM(alpha,A,B,C);
      }
    }
#ifdef XDEBUG
    if (Norm(C-C2) > 0.001*max(RealType(T)(1),Norm(C)+Norm(A)+Norm(B))) {
      cerr<<"MultMM alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      cerr<<"Norm(C2-C) = "<<Norm(C2-C)<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_TriMatrixArith_F.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


