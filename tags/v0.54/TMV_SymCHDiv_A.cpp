
#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Sym.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define CH_BLOCKSIZE TMV_BLOCKSIZE
#define CH_BLOCKSIZE2 (TMV_BLOCKSIZE/2)
#else
#define CH_BLOCKSIZE 64
#define CH_BLOCKSIZE2 2
#endif

  //
  // Decompose
  //

  template <bool cm, class T> inline void NonBlockHermCH_Decompose(
      const SymMatrixView<T>& A)
  {
    // Cholesky Decompostion is basically a simpler version of LU Decomp,
    // since U is just Lt for Hermitian matrices.  Cholesky cannot
    // be used for complex symmetric matrices.
    //
    // Also, the matrix must be positive definite.
    // Positive definite by definition means that:
    // - All eignevalues are positive.
    // Necessary, but not sufficient, conditions of positive definite
    // matrices are:
    // - All diagonal elements are positive.
    // - Aii + Ajj > 2*Re[Aij]
    // - Det(A) > 0
    // - MaxElement(A) is on main diagonal
    // A necessary and sufficient condition is:
    // - det(A(0:i,0:i)) > 0 for all 0<i<=N
    //
    // We want to decompose the matrix (input as A) into L * Lt.
    //
    // The equations for the rows of A are:
    //
    // A(i,0:i)   = L(i,0:i) * Lt(0:i,0:i)
    // A(i,i)     = L(i,0:i) * Lt(0:i,i) + |L(i,i)|^2
    // A(i,i+1:N) = L(i,0:i) * Lt(0:i,i+1:N) + L(i,i) * Lt(i,i+1:N)
    // 
    // These are solved as: (Taking A to be stored in its LowerTriangle part)
    //
    // L(i,0:i) = A(i,0:i) % L(0:i,0:i)t
    // L(i,i) = sqrt( A(i,i) - NormSq(L(i,0:i)) )
    // [ ie. row by row from the top. ]
    //
    // or 
    // 
    // L(i,i) = sqrt( A(i,i) - NormSq(L(i,0:i)) )
    // L(i+1:N,i) = (A(i+1:N,i) - L(i,0:i) * L(i+1:N,0:i)t) / L(i,i)
    // [ ie. col by col from the left ]
    //
    // The second version is not so good as written, since it accesses L
    // by both rows and columns, which means there will be non-unit strided
    // vectors in the algorithm.  However, we can save it by
    // noting that the second equation can be calculated piecemeal:
    //
    // L(i+1:N,i) = A(i+1:N,i)
    // for(j=0..i) L(i+1:N,i) -= L(i,j) * L(i+1:N,j)*
    // L(i+1:N,i) /= L(i,i)
    //
    // (The row access in the first equation can proceed similarly.)
    // 
    // Then note that for a given j, all the i's with i>j can have 
    // L(i+1:N,i) calculated at once.  So the final column algorithm
    // is:
    //
    // L = A.LowerTri();
    // for(j=0..N) 
    //   L(j,j) = sqrt(L(j,j))
    //   L(j+1:N,j) /= L(j,j)
    //   L(j+1:N,j+1:N) -= L(j+1:N,j) ^ L(j+1:N,j)*
    //
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(IsReal(T()) || A.isherm());
    TMVAssert(cm == A.iscm());
    const size_t N = A.size();
#ifdef XDEBUG
    Matrix<T> A0 = A;
    /*
    for(size_t i=1;i<=N;i++) {
      T d = Matrix<T>(A.SubSymMatrix(0,i)).Det();
      if (REAL(d) < 0) {
	cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
	cerr<<"Det(0.."<<i<<") = "<<d<<endl;
      }
    }
    */
#endif

    const VectorView<RealType(T)> Adiag = A.Real().diag();
    if (cm) {
      RealType(T)* Ajj= Adiag.ptr();
      const size_t ds = Adiag.step();
      for(size_t j=0;j<N-1;++j,Ajj+=ds) {
	if (*Ajj <= RealType(T)(0)) 
	  throw NonPosDefHermMatrix<T>(A);
	*Ajj = sqrt(*Ajj);
	A.col(j,j+1,N) /= *Ajj;
	A.SubSymMatrix(j+1,N) -= A.col(j,j+1,N) ^ A.col(j,j+1,N).Conjugate();
      }
      if (*Ajj <= RealType(T)(0)) 
	throw NonPosDefHermMatrix<T>(A);
      *Ajj = sqrt(*Ajj);
    } else {
      RealType(T)* Aii= Adiag.ptr();
      const size_t ds = Adiag.step();
      if (*Aii <= RealType(T)(0)) 
	throw NonPosDefHermMatrix<T>(A);
      *Aii = sqrt(*Aii);
      Aii+=ds;
      for(size_t i=1; i<N; ++i,Aii+=ds) {
	A.row(i,0,i) %= A.LowerTri().SubTriMatrix(0,i).Adjoint();
	*Aii -= NormSq(A.row(i,0,i));
	if (*Aii <= RealType(T)(0)) 
	  throw NonPosDefHermMatrix<T>(A);
	*Aii = sqrt(*Aii);
      }
    }
#ifdef XDEBUG
    RealType(T) norml = Norm(A.LowerTri());
    Matrix<T> A2 = A.LowerTri()*A.LowerTri().Adjoint();
    if (Norm(A2-A0) > 0.001*SQR(norml)) {
      cerr<<"HermCHDecomp: A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"Done: A = "<<A<<endl;
      cerr<<"L = "<<A.LowerTri()<<endl;
      cerr<<"Lt = "<<A.LowerTri().Adjoint()<<endl;
      cerr<<"L*Lt = "<<A2<<endl;
      cerr<<"norm(diff) = "<<Norm(A2-A0)<<"  Norm(L) = "<<norml<<endl;
      abort();
    }
#endif
  }

  template <bool cm, class T> inline void BlockHermCH_Decompose(
      const SymMatrixView<T>& A)
  {
    // If A is large, we can take advantage of Blas Level 3 speed
    // by partitioning the matrix into blocks.
    //
    // A = [ A00  A01 ]
    //     [ A10  A11 ]
    // A00, A11 are Hermitian,
    // A01 = A10t
    // A00 is k by k, A11 is (N-k) by (N-k)
    //
    // Then, the L Lt decomposition of A is:
    //
    // A = [ L00   0  ] [ L00t L10t ]
    //     [ L10  L11 ] [  0   L11t ]
    //   = [ L00 L00t         L00 L10t      ]
    //     [ L10 L00t   L10 L10t + L11 L11t ]
    //
    // So A00 = L00 L00t
    //    A10 = L10 L00t
    //    A11 = L10 L10t + L11 L11t
    //
    // The block routine is:
    //
    // 1) First Cholesky decompose A00 into L00 L00t.  
    // 2) Solve for L10 = A10 L00t^-1.
    // 3) Find A11' = A11 - L10 L10t.
    // 4) Repeat with A11'.
    //
#ifdef XDEBUG
    //cerr<<"Start Block HermCH_Decomp: A = "<<Type(A)<<"  "<<A<<endl;
    HermMatrix<T> A0 = A;
    HermMatrix<T,Lower,ColMajor> A2 = A;
    NonBlockHermCH_Decompose(A2.View());
#endif
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(IsReal(T()) || A.isherm());
    TMVAssert(cm == A.iscm());

    const size_t N = A.size();

    for (size_t jk=0; jk<N; jk+=CH_BLOCKSIZE)
    {
      size_t jkpk = min(jk+CH_BLOCKSIZE,N);

      SymMatrixView<T> A00 = A.SubSymMatrix(jk,jkpk);
      SymMatrixView<T> A11 = A.SubSymMatrix(jkpk,N);
      MatrixView<T> A10 = A.SubMatrix(jkpk,N,jk,jkpk);

      NonBlockHermCH_Decompose<cm>(A00);

      UpperTriMatrixView<T> L00t = A00.UpperTri();
      A10 %= L00t;

      A11 -= A10 * A10.Adjoint();
    }

#ifdef XDEBUG
    ConstLowerTriMatrixView<T> L = A.LowerTri();
    Matrix<T> AA = L * L.Adjoint();
    if (Norm(AA-A0) > 0.001*Norm(A0)) {
      cerr<<"Done Block CH: \n";
      cerr<<"A (block) = "<<A<<endl;
      cerr<<"A (nonblock) = "<<A2<<endl;
      cerr<<"Orig A = "<<A0<<endl;
      cerr<<"LLt (block) = "<<AA<<endl;
      cerr<<"Norm(A-A2) = "<<Norm(A-A2)<<endl;
      cerr<<"Norm(L-L2) = "<<Norm(L-A2.LowerTri())<<endl;
      cerr<<"Norm(AA-A0) = "<<Norm(AA-A0)<<endl;
      abort();
    }
#endif
  }

  template <bool cm, class T> inline void RecursiveHermCH_Decompose(
      const SymMatrixView<T>& A)
  {
#ifdef XDEBUG
    //cerr<<"Start Recursive HermCH_Decomp: A = "<<Type(A)<<"  "<<A<<endl;
    HermMatrix<T> A0 = A;
    HermMatrix<T,Lower,ColMajor> A2 = A;
    NonBlockHermCH_Decompose<true>(A2.View());
#endif
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.isherm());
    TMVAssert(cm == A.iscm());
    const size_t N = A.size();

    if (N > CH_BLOCKSIZE2) {
      size_t k = N/2;
      if (k > CH_BLOCKSIZE) k = k/CH_BLOCKSIZE*CH_BLOCKSIZE;

      SymMatrixView<T> A00 = A.SubSymMatrix(0,k);
      SymMatrixView<T> A11 = A.SubSymMatrix(k,N);
      MatrixView<T> A10 = A.SubMatrix(k,N,0,k);

      RecursiveHermCH_Decompose<cm>(A00);

      UpperTriMatrixView<T> L00t = A00.UpperTri();
      A10 %= L00t;
      A11 -= A10 * A10.Adjoint();

      RecursiveHermCH_Decompose<cm>(A11);

    } else if (CH_BLOCKSIZE2 > 2 && N > 2) {
      NonBlockHermCH_Decompose<cm>(A);
    } else if (N > 0) {
      T* Aptr = A.ptr();
      RealType(T) A00 = REAL(*Aptr);
      if (A00 <= RealType(T)(0)) 
	throw NonPosDefHermMatrix<T>(A);
      *Aptr = SQRT(A00);

      if (N == 2) {
	A00 = REAL(*Aptr);
	if (cm) {
	  T A10 = (*(Aptr+1) /= A00);
	  T* A11ptr = Aptr+1+A.stepj();
	  RealType(T) A11 = REAL(*A11ptr);
	  A11 -= NORM(A10);
	  if (A11 <= RealType(T)(0)) 
	    throw NonPosDefHermMatrix<T>(A);
	  *(A11ptr) = SQRT(A11);
	} else {
	  const int si = A.stepi();
	  T A10 = (*(Aptr+si) /= A00);
	  T* A11ptr = Aptr+1+si;
	  RealType(T) A11 = REAL(*A11ptr);
	  A11 -= NORM(A10);
	  if (A11 <= RealType(T)(0)) 
	    throw NonPosDefHermMatrix<T>(A);
	  *(A11ptr) = SQRT(A11);
	}
      }
    }
#ifdef XDEBUG
    ConstLowerTriMatrixView<T> L = A.LowerTri();
    Matrix<T> AA = L * L.Adjoint();
    if (Norm(AA-A0) > 0.001*Norm(A0)) {
      cerr<<"Done Recursive CH: \n";
      cerr<<"A (block) = "<<A<<endl;
      cerr<<"A (nonblock) = "<<A2<<endl;
      cerr<<"Orig A = "<<A0<<endl;
      cerr<<"LLt (block) = "<<AA<<endl;
      cerr<<"Norm(A-A2) = "<<Norm(A-A2)<<endl;
      cerr<<"Norm(L-L2) = "<<Norm(L-A2.LowerTri())<<endl;
      cerr<<"Norm(AA-A0) = "<<Norm(AA-A0)<<endl;
      abort();
    }
#endif
  }

  template <class T> inline void NonLapHermCH_Decompose(const SymMatrixView<T>& A)
  {
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(IsReal(T()) || A.isherm());

    if (A.iscm())
      RecursiveHermCH_Decompose<true>(A);
    else
      RecursiveHermCH_Decompose<false>(A);
    
    /*
    if (A.size() >= 2*CH_BLOCKSIZE) 
      if (A.iscm())
	BlockHermCH_Decompose<true>(A);
      else
	BlockHermCH_Decompose<false>(A);
    else
      if (A.iscm())
	NonBlockHermCH_Decompose<true>(A);
      else
	NonBlockHermCH_Decompose<false>(A);
	*/
  }

#ifdef ALAP
  template <class T> inline void LapHermCH_Decompose(
      const SymMatrixView<T>& A)
  { NonLapHermCH_Decompose(A); }
  template <> inline void LapHermCH_Decompose(const SymMatrixView<double>& A)
  {
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(A.ct()==NonConj);

    int n = A.size();
    int lda = A.iscm() ? A.stepj() : A.stepi();

    LAPNAME(dpotrf) (LAPCM A.iscm()?LAPCH_LO:LAPCH_UP, LAPV(n),
	LAPP(A.ptr()),LAPV(lda) LAPINFO LAP1);
    if (Lap_info > 0) throw NonPosDefHermMatrix<double>(A);
    LAP_Results("dpotrf");
  }
  template <> inline void LapHermCH_Decompose(
      const SymMatrixView<complex<double> >& A)
  {
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.isherm());

    int n = A.size();
    int lda = A.iscm() ? A.stepj() : A.stepi();
    LAPNAME(zpotrf) (LAPCM A.iscm()?LAPCH_LO:LAPCH_UP, LAPV(n),
	LAPP(A.ptr()),LAPV(lda) LAPINFO LAP1);
    if (Lap_info > 0) throw NonPosDefHermMatrix<complex<double> >(A);
    LAP_Results("zpotrf");
  }
#ifndef NOFLOAT
  template <> inline void LapHermCH_Decompose(const SymMatrixView<float>& A)
  {
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(A.ct()==NonConj);

    int n = A.size();
    int lda = A.iscm() ? A.stepj() : A.stepi();
    LAPNAME(spotrf) (LAPCM A.iscm()?LAPCH_LO:LAPCH_UP, LAPV(n),
	LAPP(A.ptr()),LAPV(lda) LAPINFO LAP1);
    if (Lap_info > 0) throw NonPosDefHermMatrix<float>(A);
    LAP_Results("spotrf");
  }
  template <> inline void LapHermCH_Decompose(
      const SymMatrixView<complex<float> >& A)
  {
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.isherm());

    int n = A.size();
    int lda = A.iscm() ? A.stepj() : A.stepi();
    LAPNAME(cpotrf) (LAPCM A.iscm()?LAPCH_LO:LAPCH_UP, LAPV(n),
	LAPP(A.ptr()),LAPV(lda) LAPINFO LAP1);
    if (Lap_info > 0) throw NonPosDefHermMatrix<complex<float> >(A);
    LAP_Results("cpotrf");
  }
#endif // NOFLOAT
#endif // ALAP
  template <class T> void HermCH_Decompose(const SymMatrixView<T>& A)
  {
    TMVAssert(IsReal(T()) || A.isherm());
    TMVAssert(A.isrm() || A.iscm());

    if (A.uplo() == Upper) HermCH_Decompose(A.Adjoint());
    else if (A.isconj()) HermCH_Decompose(A.Conjugate());
    else if (A.size() > 0) {
#ifdef ALAP 
      LapHermCH_Decompose(A);
#else
      NonLapHermCH_Decompose(A);
#endif
    }
  }

#define InstFile "TMV_SymCHDiv_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


