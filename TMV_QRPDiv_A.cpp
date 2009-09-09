
#include "TMV.h"
#include "TMV_Householder.h"
#include "TMV_Tri.h"
#include "TMV_VectorAux.h"

//#define XDEBUG

namespace tmv {

  bool StrictQRP = false;

#ifdef TMV_BLOCKSIZE
#define QRP_BLOCKSIZE TMV_BLOCKSIZE
#else
#define QRP_BLOCKSIZE 32
#endif

  //
  // QRP Decompose
  //

  template <class T> inline void NonBlockQRP_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& beta, size_t* P, T& det)
  {
    // Decompose A into A = Q R P
    // where Q is unitary, R is upper triangular, and P is a permutation
    // Q and R are stored in the same matrix (output of A), 
    // with the beta's for the Householder matrices returned in beta.
    
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(beta.size() == A.rowsize());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);
    const size_t M = A.colsize();
    const size_t N = A.rowsize();
    const int Astepj = A.stepj();

    // Keep track of the norm of each column
    // When considering column j, these are actually just the norm
    // of each column from j:M, not 0:M.
    Vector<RealType(T)> colnormsq(N);
    for(size_t j=0;j<N;++j) colnormsq(j) = NormSq(A.col(j));
    RealType(T) anormsq = colnormsq.SumElements();
    RealType(T) thresh = RealType(T)(N) * SQR(Epsilon<T>()) * anormsq;
    //cerr<<"thresh = "<<thresh<<endl;
    // Set to 0 any diag element whose norm is < epsilon * |A|
    RealType(T) recalcthresh(0);
    // recalcthresh is the threshold for recalculating the norm to account for
    // rounding errors in the subtractions which keep track of it..
    // The is set to sqrt(Epsilon) * the largest normsq whenever we 
    // recalculate the norms. 

    for(size_t j=0;j<N;++j) {
      //cerr<<"j = "<<j<<" colnormsq = "<<colnormsq(j)<<endl;

      if (StrictQRP || j==0 || colnormsq(j) < recalcthresh) {
	// Find the column with the largest norm
	size_t jpiv;
	RealType(T) maxnormsq = colnormsq.SubVector(j,N).MaxElement(&jpiv);
	if (j==0) recalcthresh = 4*SqrtEpsilon<T>() * maxnormsq;
	// Note: jpiv is relative to the SubVector(j,N)

	//cerr<<"jpiv = "<<jpiv+j<<", maxnormsq = "<<maxnormsq<<endl;

	// If the largest colnormsq is lower than the recalulation threshold,
	// then recalc all colnormsq's, and redetermine max.
	if (maxnormsq < recalcthresh) {
	  for(size_t k=j;k<N;++k) colnormsq(k) = NormSq(A.col(k,j,M));
	  maxnormsq = colnormsq.SubVector(j,N).MaxElement(&jpiv);
	  recalcthresh = 4*SqrtEpsilon<T>() * maxnormsq;
	  if (recalcthresh < thresh) recalcthresh = thresh;
	}

	// If maxnormsq = 0 (technically < thresh to account for rounding)
	// then the rest of the R matrix is 0, and the Householder matrices 
	// are identities (indicated by 0's in the Q part of the matrix).
	if (maxnormsq < thresh) {
	  //cerr<<"Do Zero\n";
	  A.SubMatrix(j,M,j,N).Zero();
	  // Already essentially zero - make it exact
	  beta.SubVector(j,N).Zero();
	  // Set the Householder matrices for these to identities
	  for(;j<N;j++) P[j] = j;
	  break;
	} else {

	  // Swap the column with the largest norm into the current column
	  if (jpiv != 0) {
	    // Add j to get real index
	    jpiv += j;
	    TMVAssert(jpiv < A.rowsize());
	    colnormsq.Swap(j,jpiv);
	    A.SwapCols(j,jpiv);
	    det = -det;
	    P[j] = jpiv;
	  } else {
	    P[j] = j;
	  }
	}
      } else P[j] = j;

      // Apply the Householder Reflection for this column
      beta(j) = Householder_Reflect(A.SubMatrix(j,M,j,N),det);

      // And update the norms for use with the next column
      const T* Ajk = A.row(j,j+1,N).cptr();
      for(size_t k=j+1;k<N;++k,Ajk+=Astepj) {
	colnormsq(k) -= NORM(*Ajk);
      }
    }
  }

  template <class T> inline void StrictBlockQRP_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& beta, size_t* P, T& det)
  {
    // Decompose A (input as A) into A = Q R P
    // where Q is unitary, R is upper triangular, and P is a permutation
    // Q and R are stored in the same matrix (A), with the beta's for
    // the Householder matrices returned in beta.
    
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(beta.size() == A.rowsize());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);

#ifdef XDEBUG
    Matrix<T> A0 = A;
#endif

    const size_t M = A.colsize();
    const size_t N = A.rowsize();
    const int Astepj = A.stepj();

    Vector<RealType(T)> colnormsq(N);
    for(size_t j=0;j<N;++j) colnormsq(j) = NormSq(A.col(j));
    RealType(T) anormsq = colnormsq.SumElements();
    RealType(T) thresh = RealType(T)(N) * SQR(Epsilon<T>()) * anormsq;
    RealType(T) recalcthresh(0);

    Matrix<T,RowMajor> ZYtA(min(size_t(QRP_BLOCKSIZE),M),N);
    // We keep track of the product ZYtA(j1:M,j1:N) [ stored as ZYtA ]
    // since this is the product that we need.  We update this one 
    // row at a time.
    
    for(size_t j1=0;j1<N;) {
      size_t j2 = min(N,j1+QRP_BLOCKSIZE);

      for(size_t j=j1,jmj1=0; j<j2; ++j,++jmj1) {
	size_t jpiv;
	RealType(T) maxnormsq = colnormsq.SubVector(j,N).MaxElement(&jpiv);
	if (recalcthresh == RealType(T)(0)) 
	  recalcthresh = 4*SqrtEpsilon<T>() * maxnormsq;

	if (maxnormsq < recalcthresh) {
	  for(size_t k=j;k<N;++k) colnormsq(k) = NormSq(A.col(k,j,M));
	  recalcthresh = RealType(T)(0);
	  j2 = j;
	} else if (maxnormsq < thresh) {
	  if (j==j1) {
	    // If first in set, then just zero the rest out and 
	    // indicate that we are done.
	    A.SubMatrix(j,M,j,N).Zero();
	    beta.SubVector(j,N).Zero();
	    for(;j<N;j++) P[j] = j;
	    j2 = N; 
	  } else {
	    // Otherwise, do the block Householder transforms on the 
	    // previous columns first.  (The next time through the block
	    // loop should still result in maxnormsq < thresh.)
	    j2 = j; 
	  }
	} else {

	  // Pivot
	  if (jpiv != 0) {
	    jpiv += j;
	    TMVAssert(jpiv < A.rowsize());
	    colnormsq.Swap(j,jpiv);
	    A.SwapCols(j,jpiv);
	    ZYtA.Rows(0,jmj1).SwapCols(j,jpiv);
	    det = -det;
	    P[j] = jpiv;
	  } else {
	    P[j] = j;
	  }

	  // Update the pivot column with Block Householder so far:
	  // A(j1:M,j) -= Y Z Yt A(j1:M,j)
	  // A(j1:j,j) has already been updated, so we only need to do A(j:M,j)
	  // A(j:M,j) -= Y(j:M,0:j) (ZYtA)(0:j,j)
	  A.col(j,j,M) -= A.SubMatrix(j,M,j1,j) * ZYtA.col(j,0,jmj1);

	  // Find Householder matrix for this column
	  beta(j) = Householder_Reflect(A.col(j,j,M),det);

	  // Update ZYtA:
	  if (beta(j) != T(0)) {
	    // (I-beta v vt)(I-Y Z Yt) A
	    // = I - Y (ZYtA) - v (beta vt A) + v (beta vt Y (ZYtA))
	    // The augmented Y now includes v in the j column, 
	    // so the augmented ZYtA now has to include in the j row:
	    // beta (vt A - vt Y ZYtA)
	    VectorView<T> vt = A.col(j,j+1,M).Conjugate();
	    // Remember, this doesn't include the implicit 1 at the top of v.
	    ZYtA.row(jmj1,j1,j+1).Zero();
	    ZYtA.row(jmj1,j+1,N) = vt * A.SubMatrix(j+1,M,j+1,N);
	    ZYtA.row(jmj1,j+1,N) += A.row(j,j+1,N);
	    Vector<T> vtY = vt * A.SubMatrix(j+1,M,j1,j);
	    vtY += A.row(j,j1,j);
	    ZYtA.row(jmj1,j1,N) -= vtY * ZYtA.SubMatrix(0,jmj1,j1,N);
	    ZYtA.row(jmj1,j1,N) *= beta(j);
	  } else ZYtA.row(jmj1,j1,N).Zero();

	  // Update row j of the rest of the matrix:
	  // A(j,j+1:N) -= (Y ZYtA)(j,j+1:N) = Y(j,j1:j+1) ZYtA(j1:j+1,j+1:N)
	  VectorView<T> Arowj = A.row(j,j+1,N);
	  Arowj -= A.row(j,j1,j)*ZYtA.SubMatrix(0,jmj1,j+1,N);
	  Arowj -= ZYtA.row(jmj1,j+1,N);

	  // Update the colnormsq values
	  const T* Ajk = Arowj.cptr();
	  for(size_t k=j+1;k<N;++k,Ajk+=Astepj) colnormsq(k) -= tmv::NORM(*Ajk);
	}
      }
      // Do the Block Householder update of the rest of the matrix:
      // A(j2:M,j2:N) -= Y(j2:M,j1:j2) ZYtA(j1:j2,j1:N)
      A.SubMatrix(j2,M,j2,N) -= A.SubMatrix(j2,M,j1,j2) * 
	ZYtA.SubMatrix(0,j2-j1,j2,N);
      j1 = j2;
    }

#ifdef XDEBUG
    Matrix<T> Q = A;
    GetQFromQR(Q.View(),beta);
    Matrix<T> AA = Q*UpperTriMatrixViewOf(A);
    AA.ReversePermuteCols(P);
    if (Norm(AA-A0) > 0.001*Norm(A0)) {
      cerr<<"StrictBlockQRP_Decompose: A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> "<<A<<endl;
      cerr<<"beta = "<<beta<<endl;
      cerr<<"P = ";
      for(size_t i=0;i<A.rowsize();i++) cerr<<P[i]<<" ";
      cerr<<endl;
      cerr<<"QRP = "<<AA<<endl;
      abort(); 
    }
#endif
  }

  template <class T> inline void MoveLowColsToEnd(
      Vector<RealType(T)>& colnormsq, RealType(T) thresh,
      size_t j1, size_t& j2, size_t& j3, const MatrixView<T>& A, size_t* P)
  {
    // Move all columns of A (whose norms are in colnormsq) with norms
    // less than thresh to the end.  j1 is the first column we need to 
    // look at.  j2 is one past the last column we need to look at.
    // j3 is one past the last good column overall.
    // On output j2 and j3 are updated lower if necessary.
    
    --j3; // temporarily, j3 is the last good column.
    while (j3 > j1 && colnormsq(j3) < thresh) --j3;
    if (j3==j1) {
      if (colnormsq(j1) < thresh) j2 = j1;
      else { j2 = j3 = j1+1; P[j1] = j1; }
      return;
    }
    if (j3 < j2) j2 = j3+1;
    for(size_t i=j1;i<j2;++i) {
      if (colnormsq(i) < thresh) {
	//cerr<<"  "<<i<<" is low ("<<colnormsq(i)<<") - swap with "<<j3<<endl;
	TMVAssert(j3 < A.rowsize());
	A.SwapCols(i,j3);
	colnormsq.Swap(i,j3);
	P[i] = j3;
	while (colnormsq(--j3) < thresh);
	if (j3 < j2) j2 = j3+1;
      } else P[i] = i;
    }
    ++j3; // j3 back to being first bad column.
  }

#ifdef XDEBUG
  inline void CheckIndex(const Vector<double>& index, const size_t* P, size_t j1)
  {
    const size_t N = index.size();
    Vector<double> index2(N);
    for(size_t k=0;k<N;k++) index2(k) = double(k);
    for(size_t k=0;k<j1;k++) index2.Swap(k,P[k]);
    if (Norm(index-index2) > 0.01) {
      cerr<<"index = "<<index<<endl;
      cerr<<"index2 = "<<index2<<endl;
      cerr<<"norm(diff) = "<<Norm(index-index2)<<endl;
      abort();
    }
  }
#endif

  template <class T> inline void LooseBlockQRP_Decompose(
      const MatrixView<T>& A, const VectorView<T>& beta, size_t* P, T& det)
  {
    // Decompose A (input as A) into A = Q R P
    // where Q is unitary, R is upper triangular, and P is a permutation
    // Q and R are stored in the same matrix (A), with the beta's for
    // the Householder matrices returned in beta.
    //
    // This loose version doesn't sort the diagonal of R exactly.
    // It only sorts them enough to make sure the 0's fall at the end.
    
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(beta.size() == A.rowsize());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);

#ifdef XDEBUG
    Matrix<T> A0 = A;
#endif

    const size_t M = A.colsize();
    const size_t N = A.rowsize();
    const int Astepj = A.stepj();

    Vector<RealType(T)> colnormsq(N);
    for(size_t j=0;j<N;++j) colnormsq(j) = NormSq(A.col(j));
    //cerr<<"colnormsq = "<<colnormsq<<endl;
    RealType(T) anormsq = colnormsq.SumElements();
    //cerr<<"anormsq = "<<anormsq<<endl;
    //cerr<<"eps = "<<Epsilon<T>()<<endl;
    RealType(T) thresh = RealType(T)(N) * SQR(Epsilon<T>()) * anormsq;
    //cerr<<"thresh = "<<thresh<<endl;

#ifdef XDEBUG
    Vector<double> index(N);
    for(size_t k=0;k<N;k++) index(k)=double(k);
#endif

    for (size_t j1 = 0; j1 < N;) {
      //cerr<<"Start outer loop: j1 = "<<j1<<endl;
      // Do as many columns as possible such that none have to have 
      // their norms recalculated.
      size_t j3=N; // j3 will be how many we have done in this loop
      // Invariant: all columns from j3..N are known to have norms that
      // need to be recalculated.
      // The recalculation is done at the end of the loop.

      size_t jpiv0;
      RealType(T) maxnormsq = colnormsq.SubVector(j1,N).MaxElement(&jpiv0);
      //cerr<<"maxnormsq = "<<maxnormsq<<endl;

      if (maxnormsq < thresh) {
	//cerr<<"Do Zero\n";
	// Zero the rest out and we are done.
	A.SubMatrix(j1,M,j1,N).Zero();
	beta.SubVector(j1,N).Zero();
	for(;j1<N;j1++) P[j1] = j1;
	break;
      } 

      // Move max column to the front:
      if (jpiv0 != 0) {
	jpiv0 += j1;
	//cerr<<"swap j1 with jpiv0 = "<<jpiv0<<endl;
	TMVAssert(jpiv0 < A.rowsize());
	A.SwapCols(j1,jpiv0);
	colnormsq.Swap(j1,jpiv0);
	P[j1] = jpiv0;
#ifdef XDEBUG
	index.Swap(j1,jpiv0);
#endif
      } else P[j1] = j1;
#ifdef XDEBUG
      CheckIndex(index,P,j1+1);
#endif

      RealType(T) recalcthresh = RealType(T)(N)*SqrtEpsilon<T>()*maxnormsq;
      if (recalcthresh < thresh) recalcthresh = thresh;
      //cerr<<"recalcthresh = "<<recalcthresh<<endl;

      TMVAssert(j1<j3);
      size_t j1x = j1+1; 
      // The first pass through, we don't want to include j1 in the 
      // MoveLowColsToEnd call.

      // Work on this one block at a time:
      while (j1 < j3) {
	//cerr<<"Start inner loop: j1 = "<<j1<<", j3 = "<<j3<<endl;
	size_t j2 = min(j3,j1+QRP_BLOCKSIZE);
	//cerr<<"j2 = "<<j2<<endl;
	TMVAssert(j1 < j2);
	MoveLowColsToEnd(colnormsq,recalcthresh,j1x,j2,j3,A,P);
	//cerr<<"After MoveLowColsToEnd: \n";
	//cerr<<"j2, j3 = "<<j2<<','<<j3<<endl;
#ifdef XDEBUG
	for(size_t k=j1x;k<j2;k++) index.Swap(k,P[k]);
	CheckIndex(index,P,j2);
#endif

	size_t origj2 = j2;
	UpperTriMatrix<T,NonUnitDiag,ColMajor> Z(j2-j1);
	//cerr<<"Block "<<j1<<".."<<j2<<endl;

	//cerr<<"colnormsq = "<<colnormsq.SubVector(j1,j2)<<endl;

	for(size_t j=j1; j<j2; ++j) {

	  if (colnormsq(j) < recalcthresh) {
#ifdef XDEBUG
	    CheckIndex(index,P,origj2);
#endif
	    --j2;
	    //cerr<<" "<<j<<" is low ("<<colnormsq(j)<<") - swap with "<<j2<<endl;
	    //cerr<<" Current Ps are "<<P[j]<<"  "<<P[j2]<<endl;
	    if (j==j2) break;
	    A.SwapCols(j,j2);
	    colnormsq.Swap(j,j2);
#ifdef XDEBUG
	    index.Swap(j,j2);
#endif
	    if (P[j2] > j2) {
	      if (P[j] > j2) {
		TMVAssert(P[j] < A.rowsize());
		TMVAssert(P[j2] < A.rowsize());
		swap(P[j],P[j2]);
		A.SwapCols(P[j],P[j2]);
		colnormsq.Swap(P[j],P[j2]);
#ifdef XDEBUG
		index.Swap(P[j],P[j2]);
#endif
	      } else {
		P[j] = P[j2];
	      }
	    } else {
	      if (P[j] > j2) P[j2] = P[j];
	      P[j] = j2;
	    }
#ifdef XDEBUG
	    CheckIndex(index,P,origj2);
#endif
	  }

	  // Find Householder matrix for this column
	  // This multiplies through to the end of the original block.
	  // This way, when we are done, the whole block has had the
	  // same Householder reflections applied to it.
	  beta(j) = Householder_Reflect(A.SubMatrix(j,M,j,origj2),det);

	  // Update Z:
	  BlockHouseholder_Augment(A.SubMatrix(j1,M,j1,j+1),
	      Z.SubTriMatrix(0,j-j1+1),CONJ(beta(j)));

	  // Update the colnormsq values within this block
	  // (No need to go all the way to origj2, since the j2..origj2
	  // columns are those with low norm already - we don't need
	  // those values until we recalculate them from scratch anyway.)
	  const T* Ajk = A.row(j,j+1,j2).cptr();
	  for(size_t k=j+1;k<j2;++k,Ajk+=Astepj) 
	    colnormsq(k) -= tmv::NORM(*Ajk);
	}
	//cerr<<"done: j2 = "<<j2<<endl;

	if (j1 < j2) {

	  // Do the Block Householder update of the rest of the matrix:
	  BlockHouseholder_LDiv(A.SubMatrix(j1,M,j1,j2),
	      Z.SubTriMatrix(0,j2-j1),A.SubMatrix(j1,M,origj2,N));

	  // Update the colnormsq values for the rest of the matrix:
	  if (M-j2 > j2-j1)
	    for(size_t k=origj2;k<N;++k) colnormsq(k) -= 
	      A.col(k,j1,j2).NormSq();
	  else 
	    for(size_t k=origj2;k<N;++k) colnormsq(k) = 
	      A.col(k,j2,M).NormSq();
	}

	if (j2 < origj2) {
	  //cerr<<"j2 = "<<j2<<", origj2 = "<<origj2<<endl;
	  //cerr<<"P["<<j2<<".."<<origj2<<"] = ";
	  //for(size_t j=j2; j<origj2; ++j) cerr<<P[j]<<" ";
	  
#ifdef XDEBUG
	  CheckIndex(index,P,origj2);
#endif
	  // Put the bad columns back where they started before this loop:
	  for(size_t j=j2; j<origj2; ++j) if (P[j] > j2) {
	    //cerr<<"Sending "<<j<<" back to "<<P[j]<<endl;
	    TMVAssert(P[j] < A.rowsize());
	    A.SwapCols(j,P[j]);
	    colnormsq.Swap(j,P[j]);
#ifdef XDEBUG
	    index.Swap(j,P[j]);
#endif
	  }
#ifdef XDEBUG
	  CheckIndex(index,P,j2);
#endif
	}

	j1 = j1x = j2;
	//cerr<<"For next loop: j1 = "<<j1<<", j2 = "<<j2<<endl;
      }
      //cerr<<"Done main loop pass\n";

      if (j3 < N) {
	// Then need to recalculate some of the colnorms:
	for(size_t k=j3;k<N;++k) colnormsq(k) = NormSq(A.col(k,j3,M));
	//cerr<<"New colnorms are "<<colnormsq.SubVector(j3,N)<<endl;
      }
    }

    if (det != T(0)) {
      for(size_t i=0;i<N;++i) if (P[i] != i) det = -det;
    }

#ifdef XDEBUG
    CheckIndex(index,P,N);
    Matrix<T> Q = A;
    GetQFromQR(Q.View(),beta);
    Matrix<T> AA = Q*UpperTriMatrixViewOf(A);
    AA.ReversePermuteCols(P);
    if (Norm(AA-A0) > 0.001*max(RealType(T)(1),Norm(A0))) {
      cerr<<"LooseBlockQRP_Decompose: \n";
      cerr<<"A = "<<Type(A)<<endl;
      if (N < 100) {
	cerr<<"  "<<A0<<endl;
	cerr<<"-> "<<A<<endl;
	cerr<<"beta = "<<beta<<endl;
	cerr<<"P = ";
	for(size_t i=0;i<N;i++) cerr<<P[i]<<" ";
	cerr<<endl;
	cerr<<"QRP = "<<AA<<endl;
	Matrix<T> diff = AA-A0;
	diff.Clip(0.0001);
	cerr<<"diff = "<<diff<<endl;
      }
      cerr<<"Rdiag = "<<A.diag()<<endl;
      cerr<<"Norm(A-QRP) = "<<Norm(AA-A0)<<endl;
      abort(); 
    }
#endif
  }

  template <class T> inline void NonLapQRP_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& beta, size_t* P, T& det)
  {
    // Decompose A (input as A) into A = Q R P
    // where Q is unitary, R is upper triangular, and P is a permutation
    // Q and R are stored in the same matrix (A), with the beta's for
    // the Householder matrices returned in beta.

    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(beta.size() == A.rowsize());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);

    if (A.rowsize() > QRP_BLOCKSIZE)
      if (StrictQRP)
	StrictBlockQRP_Decompose(A,beta,P,det);
      else
	LooseBlockQRP_Decompose(A,beta,P,det);
    else
      NonBlockQRP_Decompose(A,beta,P,det);
  }

#ifdef LAP
  template <class T> inline void LapQRP_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& beta, size_t* P, T& det)
  { NonLapQRP_Decompose(A,beta,P,det); }
//#ifndef MKL
  // With the Intel MKL version, neither LAP QRP routine (the other,
  // non-level-3 version is geqpf) seems to do any actual pivoting.
  // So until I figure out why and how to get them to work, just
  // do the NonLap routine.
  template <> inline void LapQRP_Decompose(const MatrixView<double>& A,
      const VectorView<double>& beta, size_t* P, double& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);
    int m = A.colsize();
    int n = A.rowsize();
    auto_array<int> lap_p(new int[n]);
    for(int i=0;i<n;++i) (lap_p.get())[i] = 0;
#ifndef LAPNOWORK
    int lwork = 3*n+1;
    double* work = LAP_DWork(lwork);
#endif
    int lda = A.stepj();
    LAPNAME(dgeqp3) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
	LAPP(lap_p.get()),LAPP(beta.ptr()) LAPWK(work) LAPVWK(lwork) LAPINFO);
    LAP_Results("dgeqp3");
    double thresh = Epsilon<double>()*A.NormF();
    for(size_t i=0;i<A.rowsize();++i) {
      if (abs(A(i,i)) < thresh) {
	TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
	A.SubMatrix(i,A.colsize(),i,A.rowsize()).Zero();
	beta.SubVector(i,beta.size()).Zero();
	break;
      }
    }
#ifndef CLAP
    for(size_t i=0;i<beta.size();++i) {
      --((lap_p.get())[i]);
    }
#endif
    ConvertIndexToPermute(n,lap_p.get(),P);
    for(size_t i=0;i<beta.size();++i) {
      if (det) {
	if (beta(i) != 0.) det = -det;
	if (P[i] != i) det = -det;
      }
    }
  }
  template <> inline void LapQRP_Decompose(
      const MatrixView<complex<double> >& A,
      const VectorView<complex<double> >& beta, size_t* P,
      complex<double>& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);
    int m = A.colsize();
    int n = A.rowsize();
    auto_array<int> lap_p(new int[n]);
    for(int i=0;i<n;++i) (lap_p.get())[i] = 0;
#ifndef LAPNOWORK
    int lwork = 2*n;
    double* rwork = LAP_DWork(lwork);
    lwork = n+1;
    complex<double>* work = LAP_ZWork(lwork);
#endif
    int lda = A.stepj();
    LAPNAME(zgeqp3) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
	LAPP(lap_p.get()),LAPP(beta.ptr()) 
	LAPWK(work) LAPVWK(lwork) LAPWK(rwork) LAPINFO);
    LAP_Results("zgeqp3");
    beta.ConjugateSelf();
    double thresh = Epsilon<double>()*A.NormF();
    for(size_t i=0;i<A.rowsize();++i) {
      TMVAssert(abs(imag(A(i,i))) < Epsilon<double>());
      if (abs(real(A(i,i))) < thresh) {
	TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
	A.SubMatrix(i,A.colsize(),i,A.rowsize()).Zero();
	beta.SubVector(i,beta.size()).Zero();
	break;
      }
    }
#ifndef CLAP
    for(size_t i=0;i<beta.size();++i) {
      --((lap_p.get())[i]);
    }
#endif
    ConvertIndexToPermute(n,lap_p.get(),P);
    for(size_t i=0;i<beta.size();++i) {
      if (det!=double(0)) {
	if (imag(beta(i)) != 0.) 
	  det *= -conj(beta(i)*beta(i))/norm(beta(i));
	else if (real(beta(i)) != 0.)
	  det = -det;
	if (P[i] != i) det = -det;
      }
    }
  }
#ifndef NOFLOAT
  template <> inline void LapQRP_Decompose(const MatrixView<float>& A,
      const VectorView<float>& beta, size_t* P, float& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);
    int m = A.colsize();
    int n = A.rowsize();
    auto_array<int> lap_p(new int[n]);
    for(int i=0;i<n;++i) (lap_p.get())[i] = 0;
#ifndef LAPNOWORK
    int lwork = 3*n+1;
    float* work = LAP_SWork(lwork);
#endif
    int lda = A.stepj();
    LAPNAME(sgeqp3) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
	LAPP(lap_p.get()),LAPP(beta.ptr()) LAPWK(work) LAPVWK(lwork) LAPINFO);
    LAP_Results("sgeqp3");
    float thresh = Epsilon<float>()*A.NormF();
    for(size_t i=0;i<A.rowsize();++i) {
      if (abs(A(i,i)) < thresh) {
	TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
	A.SubMatrix(i,A.colsize(),i,A.rowsize()).Zero();
	beta.SubVector(i,beta.size()).Zero();
	break;
      }
    }
#ifndef CLAP
    for(size_t i=0;i<beta.size();++i) {
      --((lap_p.get())[i]);
    }
#endif
    ConvertIndexToPermute(n,lap_p.get(),P);
    for(size_t i=0;i<beta.size();++i) {
      if (det) {
	if (beta(i) != 0.) det = -det;
	if (P[i] != i) det = -det;
      }
    }
  }
  template <> inline void LapQRP_Decompose(
      const MatrixView<complex<float> >& A,
      const VectorView<complex<float> >& beta, size_t* P,
      complex<float>& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.iscm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);
    int m = A.colsize();
    int n = A.rowsize();
    auto_array<int> lap_p(new int[n]);
    for(int i=0;i<n;++i) (lap_p.get())[i] = 0;
#ifndef LAPNOWORK
    int lwork = 2*n;
    float* rwork = LAP_SWork(lwork);
    lwork = n+1;
    complex<float>* work = LAP_CWork(lwork);
#endif
    int lda = A.stepj();
    LAPNAME(cgeqp3) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
	LAPP(lap_p.get()),LAPP(beta.ptr()) 
	LAPWK(work) LAPVWK(lwork) LAPWK(rwork) LAPINFO);
    LAP_Results("cgeqp3");
    float thresh = Epsilon<float>()*A.NormF();
    for(size_t i=0;i<A.rowsize();++i) {
      TMVAssert(abs(imag(A(i,i))) < Epsilon<float>());
      if (abs(real(A(i,i))) < thresh) {
	TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
	A.SubMatrix(i,A.colsize(),i,A.rowsize()).Zero();
	beta.SubVector(i,beta.size()).Zero();
	break;
      }
    }
    beta.ConjugateSelf();
#ifndef CLAP
    for(size_t i=0;i<beta.size();++i) {
      --((lap_p.get())[i]);
    }
#endif
    ConvertIndexToPermute(n,lap_p.get(),P);
    for(size_t i=0;i<beta.size();++i) {
      if (det!=float(0)) {
	if (imag(beta(i)) != 0.) 
	  det *= -conj(beta(i)*beta(i))/norm(beta(i));
	else if (real(beta(i)) != 0.)
	  det = -det;
	if (P[i] != i) det = -det;
      }
    }
  }
#endif // FLOAT
//#endif // !MKL
#endif // LAP

  template <class T> void QRP_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& beta, size_t* P, T& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);

#ifdef XDEBUG
    Matrix<T> A0 = A;
#ifdef LAP
    Matrix<T> A2 = A;
    Vector<T> beta2 = beta;
    auto_array<size_t> P2(new size_t[beta.size()]);
    T det2=det;
    NonLapQRP_Decompose(A2.View(),beta2.View(),P2.get(),det2);
#endif
#endif

    if (A.rowsize() > 0) {
#ifdef LAP
      if (A.iscm()) {
	LapQRP_Decompose(A,beta,P,det);
      } else {
	NonLapQRP_Decompose(A,beta,P,det);
      }
#else
      NonLapQRP_Decompose(A,beta,P,det);
#endif
    }
#ifdef XDEBUG
    Matrix<T> Q = A;
    GetQFromQR(Q.View(),beta);
    Matrix<T> R = UpperTriMatrixViewOf(A);
    Matrix<T> AA = Q*R;
    AA.ReversePermuteCols(P);
    if (Norm(AA-A0) > 0.001*Norm(A0)) {
      cerr<<"BlockQRP_Decompose: A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> "<<A<<endl;
      cerr<<"beta = "<<beta<<endl;
      cerr<<"P = ";
      for(size_t i=0;i<A.rowsize();i++) cerr<<P[i]<<" ";
      cerr<<endl;
      cerr<<"QRP = "<<AA<<endl;
#ifdef LAP
      cerr<<"NonLap: A = "<<A2<<endl;
      cerr<<"beta2 = "<<beta2<<endl;
      cerr<<"det2 = "<<det2<<endl;
#endif
      cerr<<"Norm(AA-A0) = "<<Norm(A0-AA)<<endl;
      abort(); 
    }
#endif
  }

  //
  // QRP Decompose - Unpacked
  //

  template <class T> void QRP_Decompose(
      const MatrixView<T>& Q, const UpperTriMatrixView<T>& R, size_t* P, T& det)
  {
    // Decompose A (input as Q) into A = Q R P
    // where Q is unitary, R is upper triangular, and P is a permutation

    const size_t N = Q.rowsize();
    TMVAssert(Q.colsize() >= N);
    TMVAssert(R.colsize() == N);
    TMVAssert(R.rowsize() == N);
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(R.ct() == NonConj);

    Vector<T> beta(N);
    QRP_Decompose(Q,beta.View(),P,det);
    R = UpperTriMatrixViewOf(Q);
    GetQFromQR(Q,beta.View());
  }

#define InstFile "TMV_QRPDiv_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


