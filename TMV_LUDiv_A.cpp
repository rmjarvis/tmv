
#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Diag.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define LU_BLOCKSIZE TMV_BLOCKSIZE
#define LU_BLOCKSIZE2 (TMV_BLOCKSIZE/2)
#else
#define LU_BLOCKSIZE 64
#define LU_BLOCKSIZE2 2
#endif

  //
  // Decompose
  //

  template <class T> inline void NonBlockLU_Decompose(
      const MatrixView<T>& A, size_t* P, T& det)
  {
    // LU Decompostion with partial pivoting.
    //
    // We want to decompose the matrix (input as A) into P * L * U
    // where P is a permutation, L is a lower triangle matrix with 1's for 
    // the diagonal, and U is an upper triangle matrix.  
    //
    // We do this calculation one column at a time.  There are other versions
    // of this algorithm which use Rank 1 updates.  This version uses mostly
    // martix-vector products and forward substitutions.  The different
    // algorithms all require the same number of operation, but this
    // one requires fewer vector touches which often means it will be 
    // slightly faster.
    //
    // After doing j cols of the calculation, we have calculated
    // U(0:j,0:j) and L(0:N,0:j) 
    // (where my a:b notation does not include the index b)
    //
    // The equation A = LU gives for the j col:
    //
    // A(0:N,j) = L(0:N,0:N) U(0:N,j)
    //
    // which breaks up into:
    //
    // (1) A(0:j,j) = L(0:j,0:N) U(0:N,j)
    // (2) A(j:N,j) = L(j:N,0:N) U(0:N,j)
    //
    // The first of these (1) simplifies to:
    // 
    // (1*) A(0:j,j) = L(0:j,0:j) U(0:j,j)
    //
    // since L is lower triangular, so L(0:j,j:N) = 0.
    // L(0:j,0:j) is already known, so this equation can be solved for
    // U(0:j,j) by forward substitution.
    // 
    // The second equation (2) simplifies to:
    //
    //      A(j:N,j) = L(j:N,0:j+1) U(0:j+1,j)
    // (2*)          = L(j:N,0:j) U(0:j,j) + L(j:N,j) U(j,j)
    //
    // since U is upper triangular so U(j+1:N,j) = 0.
    // Since we now know U(0:j,j) from (1*) above, this equation can
    // be solved for the product L(j:N,j) U(j,j)
    // 
    // This means we have some leeway on the values for L(j,j) and U(j,j),
    // as only their product is specified.
    //
    // If we take U to have unit diagonal, then L(j,j) is set here along
    // with the rest of the L(j:N,j) column.  However, this will mean that the
    // forward substutions in the (1*) steps will require divisions by the 
    // non-unit-diagonal elements of L.  It is faster to take L to have
    // unit-diagonal elements, and do the division by U(j,j) here, since then
    // we can calculate 1/U(j,j) and multiply.  So 1 division and N-j
    // multiplies which is generally faster than N-j divisions.
    //
    // However, another potential problem is that U(j,j) could be 0, or close
    // to 0.  This would lead to either an error or inaccurate results.
    // Thus we add a step in the middle of the (2*) calculation:
    //
    // Define v(j:N) = A(j:N,j) - L(j:N,0:j) U(0:j,j)
    // 
    // We search v for the element with the largest absolute value and apply
    // a permutation to swap it into the j spot.  This element then becomes 
    // U(j,j), which is then the divisor for the rest of the vector.  This 
    // will minimize the possibility of roundoff errors due to small U(j,j)'s.
    // 
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.iscm());
    const size_t N = A.rowsize();
    const size_t M = A.colsize();
    const size_t R = min(N,M);
#ifdef XDEBUG
    Matrix<T> A0 = A;
#endif

    const T* Ujj = A.cptr();
    const int Ads = A.stepj()+1;
    size_t* Pj = P;

    for (size_t j=0; j<R; ++j,Ujj+=Ads,++Pj)
    {
      if (j > 0) {
	// Solve for U(0:j,j))
	A.col(j,0,j) /= LowerTriMatrixViewOf(A.SubMatrix(0,j,0,j),UnitDiag);

	// Solve for v = L(j:M,j) U(j,j)
	A.col(j,j,M) -= A.SubMatrix(j,M,0,j) * A.col(j,0,j);
      }

      // Find the pivot element
      size_t ip;
      A.col(j,j,M).MaxAbsElement(&ip);
      // ip is relative to j index, not absolute.

      // Swap the pivot row with j if necessary
      if (ip != 0) {
	ip += j;
	TMVAssert(ip < A.colsize());
	TMVAssert(j < A.colsize());
	A.SwapRows(ip,j);  // This does both Lkb and A'
	*Pj = ip;
	det = -det;
      } else *Pj = j;

      // Solve for L(j+1:M,j)
      // If Ujj is 0, then all of the L's are 0.
      // ie. Ujj Lij = 0 for all i>j
      // Any value for Lij is valid, so leave them 0.
      if (*Ujj != T(0)) 
	A.col(j,j+1,M) /= *Ujj;
    }
    if (N > M) {
      // Solve for U(0:M,M:N))
      A.Cols(M,N) /= LowerTriMatrixViewOf(A.Cols(0,M),UnitDiag);
    }
#ifdef XDEBUG
    Matrix<T> L = A;
    UpperTriMatrixViewOf(L).OffDiag().Zero();
    L.diag().SetAllTo(T(1));
    Matrix<T> U = UpperTriMatrixViewOf(A,NonUnitDiag);
    Matrix<T> AA = L * U;
    AA.ReversePermuteRows(P,0,R);
    if (Norm(AA-A0) > 0.001*Norm(A0)) {
      cerr<<"Done NonBlock LU: \n";
      cerr<<"A0 = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"LU = "<<A<<endl;
      cerr<<"P = (";
      for(size_t i=0;i<R;i++) cerr<<P[i]<<" ";
      cerr<<")\n";
      cerr<<"AA = "<<AA<<endl;
      cerr<<"Norm(A0-AA) = "<<Norm(AA-A0)<<endl;
      abort();
    }
#endif
  }

  /*
  template <class T> inline void BlockLU_Decompose(
      const MatrixView<T>& A, size_t* P, T& det)
  {
    // If A is large, we can take advantage of Blas Level 3 speed
    // by partitioning the matrix by columns.
    //
    // We do this calculation one block at a time down the diagonal. 
    // Each block contains k columns (expect for possibly the last 
    // block which may be fewer).
    //
    // For the first block, we decompose A into:
    //
    // ( A00 A01 ) = ( L00  0  ) ( U00 U01 )
    // ( A10 A11 )   ( L10  A~ ) (  0   1  )
    //
    // From this we obtain:
    //
    // (1) A00 = L00 U00
    // (2) A10 = L10 U00
    // (3) A01 = L00 U01
    // (4) A11 = L10 U01 + A~
    //
    // For (1) we decompose A00 in place using the non-blocked algorithm.
    // (2) and (3) then give us L10 and U01.
    // Finally, (4) lets us solve for A~.
    // Repeat until done.
    //
    // With pivoting, the only real change is to combine equations (1),(2)
    // and solve both together with the non-blocked algorithm.
    //
#ifdef XDEBUG
    Matrix<T> A0 = A;
    Matrix<T,ColMajor> A2 = A;
    int P2[A.colsize()];
    T det2=1;
    NonBlockLU_Decompose(A2.View(),P2,det2);
#endif

    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.iscm());
    const size_t N = A.rowsize();
    const size_t M = A.colsize();
    const size_t R = min(N,M);

    for (size_t jk=0; jk<R; jk+=LU_BLOCKSIZE)
    {
      size_t jkpk = min(jk+LU_BLOCKSIZE,R);
      // Solve for L00, U00

      NonBlockLU_Decompose(A.SubMatrix(jk,M,jk,jkpk),P+jk,det);

      // Apply the permutation to the rest of the matrix
      if (jk > 0) {
	A.SubMatrix(jk,M,0,jk).PermuteRows(P+jk,0,jkpk-jk);
      }
      if (jkpk < N) {
	A.SubMatrix(jk,M,jkpk,N).PermuteRows(P+jk,0,jkpk-jk);

	// Solve for U01
	A.SubMatrix(jk,jkpk,jkpk,N) /= 
	  LowerTriMatrixViewOf(A.SubMatrix(jk,jkpk,jk,jkpk),UnitDiag);

	// Solve for A~
	if (jkpk < M)
	  A.SubMatrix(jkpk,M,jkpk,N) -= A.SubMatrix(jkpk,M,jk,jkpk) *
	    A.SubMatrix(jk,jkpk,jkpk,N);
      }
      for(size_t i=jk;i<jkpk;++i) P[i]+=jk;
    }
#ifdef XDEBUG
    Matrix<T> L = A;
    UpperTriMatrixViewOf(L).OffDiag().Zero();
    L.diag().SetAllTo(T(1));
    Matrix<T> U = UpperTriMatrixViewOf(A,NonUnitDiag);
    Matrix<T> AA = L * U;
    AA.ReversePermuteRows(P,0,R);
    if (Norm(AA-A0) > 0.001*Norm(A0)) {
      cerr<<"Done Block LU: \n";
      cerr<<"A0 = "<<A0<<endl;
      cerr<<"LU = "<<A<<endl;
      cerr<<"P = (";
      for(size_t i=0;i<R;i++) cerr<<P[i]<<" ";
      cerr<<")\n";
      cerr<<"A2 = "<<A2<<endl;
      cerr<<"Norm(A-A2) = "<<Norm(A-A2)<<endl;
      cerr<<"Norm(AA-A0) = "<<Norm(AA-A0)<<endl;
      cerr<<"correct P = (";
      for(size_t i=0;i<R;i++) cerr<<P2[i]<<" ";
      cerr<<")\n";
      abort();
    }
#endif
  }
*/

  template <class T> inline void RecursiveLU_Decompose(
      const MatrixView<T>& A, size_t* P, T& det)
  {
    //cerr<<"Start recursive LU: A = "<<A<<endl;
    // The recursive LU algorithm is similar to the block algorithm, except 
    // that the block is roughly half the size of the whole matrix.
    // We keep dividing the matrix in half (column-wise) until we get down
    // to an Mx2 or Mx1 matrix.

#ifdef XDEBUG
    Matrix<T> A0 = A;
    Matrix<T,ColMajor> A2 = A;
    size_t P2[A.colsize()];
    T det2=1;
    NonBlockLU_Decompose(A2.View(),P2,det2);
#endif

    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.iscm());
    const size_t N = A.rowsize();
    const size_t M = A.colsize();
    const size_t R = min(N,M);

    if (R > LU_BLOCKSIZE2) {
      //cerr<<"large R\n";
      // Split N in half, with N1 being rounded to multiple of BLOCKSIZE
      // if appropriate.
      size_t N1 = R/2;
      if (N1 > LU_BLOCKSIZE) N1 = (N1/LU_BLOCKSIZE)*LU_BLOCKSIZE;

      MatrixView<T> A0 = A.Cols(0,N1);
      MatrixView<T> A00 = A0.Rows(0,N1);
      MatrixView<T> A10 = A0.Rows(N1,M);
      MatrixView<T> A1 = A.Cols(N1,N);
      MatrixView<T> A01 = A1.Rows(0,N1);
      MatrixView<T> A11 = A1.Rows(N1,M);

      // Decompose left half into PLU
      RecursiveLU_Decompose(A0,P,det);

      // Apply the permutation to the right half of the matrix
      //cerr<<"Before A1 Permuterows\n";
      //cerr<<"N1 = "<<N1<<endl;
      //cerr<<"P(0..N1) = ";
      //for(size_t i=0;i<N1;++i) cerr<<P[i]<<" ";
      //cerr<<endl;
      A1.PermuteRows(P,0,N1);
      //cerr<<"After Permuterows\n";

      // Solve for U01
      A01 /= LowerTriMatrixViewOf(A00,UnitDiag);

      // Solve for A~
      A11 -= A10 * A01;

      // Decompose A~ into PLU
      RecursiveLU_Decompose(A11,P+N1,det);
      for(size_t i=N1;i<R;++i) P[i]+=N1;

      // Apply the new permutations to the left half
      //cerr<<"Before A0 Permuterows\n";
      //cerr<<"P(N1..R) = ";
      //for(size_t i=N1;i<R;++i) cerr<<P[i]<<" ";
      //cerr<<endl;
      A0.PermuteRows(P,N1,R);
      //cerr<<"After Permuterows\n";
    } else if (LU_BLOCKSIZE2 > 2 && R > 2) {
      NonBlockLU_Decompose(A,P,det);
    } else if (R == 2) {
      // Same as NonBlock version, but with R==2 hard coded
      VectorView<T> A0 = A.col(0);
      VectorView<T> A1 = A.col(1);

      size_t ip0,ip1;
      RealType(T) piv = A0.MaxAbsElement(&ip0);
      if (piv != RealType(T)(0)) {
	if (ip0 != 0) {
	  A0.Swap(ip0,0);
	  A1.Swap(ip0,0);
	  det = -det;
	} 

	// A0.SubVector(1,M) /= A00;
	// A1.SubVector(1,M) -= A0.SubVector(1,M) * A01;
	const T invA00 = RealType(T)(1)/(*A0.cptr());
	const T A01 = (*A1.cptr());
	piv = RealType(T)(0); // next pivot element
	ip1 = 1;
	T* Ai0 = A0.ptr()+1;
	T* Ai1 = A1.ptr()+1;
	for(size_t i=1;i<M;++i,++Ai0,++Ai1) {
	  *Ai0 *= invA00;
	  *Ai1 -= *Ai0 * A01;
	  RealType(T) absAi1 = abs(*Ai1);
	  if (absAi1 > piv) { piv = absAi1; ip1=i; }
	}
      } else {
	piv = A1.SubVector(1,M).MaxAbsElement(&ip1); 
	ip1++;
      }

      if (piv != RealType(T)(0) && M>2) {
	if (ip1 != 1) {
	  A1.Swap(ip1,1);
	  A0.Swap(ip1,1);
	  det = -det;
	} 

	const T A11 = (*(A1.cptr()+1));
	A1.SubVector(2,M) /= A11;
      }

      if (N > 2) {
	// M=2, N>2, so solve for U(0:2,2:N))
	// A.Cols(2,N).PermuteRows(P);
	if (*P == 1) A.Cols(2,N).SwapRows(0,1);
	// A.Cols(2,N) /= LowerTriMatrixViewOf(A.Cols(0,2),UnitDiag);
	const T A01 = (*A1.cptr());
	A.row(1,2,N) -= A01 * A.row(0,2,N);
      }
      P[0] = ip0;
      P[1] = ip1;
    } else if (R == 1) {
      // Same as NonBlock version, but with R==1 hard coded
      VectorView<T> A0 = A.col(0);

      RealType(T) piv = A0.MaxAbsElement(P);

      if (piv != RealType(T)(0)) {
	if (*P != 0) {
	  A0.Swap(*P,0);
	  det = -det;
	}
	A0.SubVector(1,M) /= (*A0.cptr());
      }
    }
#ifdef XDEBUG
    Matrix<T> L = A;
    UpperTriMatrixViewOf(L).OffDiag().Zero();
    L.diag().SetAllTo(T(1));
    Matrix<T> U = UpperTriMatrixViewOf(A,NonUnitDiag);
    Matrix<T> AA = L * U;
    AA.ReversePermuteRows(P,0,R);
    if (Norm(AA-A0) > 0.001*Norm(A0)) {
      cerr<<"Done Recursive LU: \n";
      cerr<<"A0 = "<<A0<<endl;
      cerr<<"LU = "<<A<<endl;
      cerr<<"P = (";
      for(size_t i=0;i<R;i++) cerr<<P[i]<<" ";
      cerr<<")\n";
      cerr<<"A2 = "<<A2<<endl;
      cerr<<"Norm(A-A2) = "<<Norm(A-A2)<<endl;
      cerr<<"Norm(AA-A0) = "<<Norm(AA-A0)<<endl;
      cerr<<"correct P = (";
      for(size_t i=0;i<R;i++) cerr<<P2[i]<<" ";
      cerr<<")\n";
      abort();
    }
#endif
  }

  template <class T> inline void NonLapLU_Decompose(
      const MatrixView<T>& A, size_t* P, T& det)
  {
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.iscm());

    RecursiveLU_Decompose(A,P,det);

    /*
    if (A.rowsize() >= 2*LU_BLOCKSIZE) 
      BlockLU_Decompose(A,P,det);
    else
      NonBlockLU_Decompose(A,P,det);
      */
  }

#ifdef ALAP
  template <class T> inline void LapLU_Decompose(
      const MatrixView<T>& A, size_t* P, T& det)
  { NonLapLU_Decompose(A,P,det); }
  template <> inline void LapLU_Decompose(
      const MatrixView<double>& A, size_t* P, double& det)
  {
    TMVAssert(A.iscm());
    TMVAssert(A.ct()==NonConj);

    int m = A.colsize();
    int n = A.rowsize();
    int lda = A.stepj();
    auto_array<int> lap_p(new int[n]);
    LAPNAME(dgetrf) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
	LAPP(lap_p.get()) LAPINFO);
    LAP_Results("dgetrf");
    for(size_t i=0;i<A.colsize();i++) {
      P[i] = (lap_p.get())[i] LAPMINUS1;
      if (P[i]!=i) det = -det;
    }
  }
  template <> inline void LapLU_Decompose(
      const MatrixView<complex<double> >& A, size_t* P,
      complex<double>& det)
  {
    TMVAssert(A.iscm());
    TMVAssert(A.ct()==NonConj);

    int m = A.colsize();
    int n = A.rowsize();
    int lda = A.stepj();
    auto_array<int> lap_p(new int[n]);
    LAPNAME(zgetrf) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
	LAPP(lap_p.get()) LAPINFO);
    LAP_Results("zgetrf");
    for(size_t i=0;i<A.colsize();i++) {
      P[i] = (lap_p.get())[i] LAPMINUS1;
      if (P[i]!=i) det = -det;
    }
  }
#ifndef NOFLOAT
#ifndef MKL
  // This is giving me a weird runtime error sometimes with MKL:
  // OMP abort: Unable to set worker thread stack size to 2098176 bytes
  // Try reducing KMP_STACKSIZE or increasing the shell stack limit.
  // So I'm cutting it out for MKL compilations
  template <> inline void LapLU_Decompose(
      const MatrixView<float>& A, size_t* P, float& det)
  {
    TMVAssert(A.iscm());
    TMVAssert(A.ct()==NonConj);

    int m = A.colsize();
    int n = A.rowsize();
    int lda = A.stepj();
    auto_array<int> lap_p(new int[n]);
    LAPNAME(sgetrf) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
	LAPP(lap_p.get()) LAPINFO);
    LAP_Results("sgetrf");
    for(size_t i=0;i<A.colsize();i++) {
      P[i] = (lap_p.get())[i] LAPMINUS1;
      if (P[i]!=i) det = -det;
    }
  }
  template <> inline void LapLU_Decompose(
      const MatrixView<complex<float> >& A, size_t* P, 
      complex<float>& det)
  {
    TMVAssert(A.iscm());
    TMVAssert(A.ct()==NonConj);

    int m = A.colsize();
    int n = A.rowsize();
    int lda = A.stepj();
    auto_array<int> lap_p(new int[n]);
    LAPNAME(cgetrf) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
	LAPP(lap_p.get()) LAPINFO);
    LAP_Results("cgetrf");
    for(size_t i=0;i<A.colsize();i++) {
      P[i] = (lap_p.get())[i] LAPMINUS1;
      if (P[i]!=i) det = -det;
    }
  }
#endif // MKL
#endif // FLOAT
#endif // ALAP

  template <class T> void LU_Decompose(
      const MatrixView<T>& A, size_t* P, T& det)
  {
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.iscm());

    if (A.colsize() > 0 && A.rowsize() > 0) {
#ifdef ALAP
      LapLU_Decompose(A,P,det);
#else
      NonLapLU_Decompose(A,P,det);
#endif
    }
  }

#define InstFile "TMV_LUDiv_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


