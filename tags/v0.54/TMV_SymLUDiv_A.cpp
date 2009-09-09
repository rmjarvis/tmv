
#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Sym.h"
#include "TMV_Diag.h"
#include "TMV_SymLUDiv_A.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define SYM_LU_BLOCKSIZE TMV_BLOCKSIZE
#else
#define SYM_LU_BLOCKSIZE 48
#endif

  //
  // Decompose
  //


  template <bool herm, class T> void SymInvert_2x2(
      T& a, T& b, T& c, T* dd)
  {
    // Invert matrix [ a  c* ] -->  1/(ab-|c|^2) [ b -c* ]
    //               [ c  b  ]                   [ -c a  ]
    if (herm) {
      RealType(T) d = REAL(a)*REAL(b) - NORM(c);
      swap(a,b);
      a/=d;
      b/=d;
      c = -c/d;
      if (dd) *dd = d;
    } else {
      T d = a*b - c*c;
      swap(a,b);
      a/=d;
      b/=d;
      c = -c/d;
      if (dd) *dd = d;
    }
  }

  // MJ: This is a bit slower for row major, even with blocking, so write the
  // row major version where A is decomposed into Lt D L rather than L D Lt.
  template <bool herm, class T> inline void NonBlockSymLU_Decompose(
      const SymMatrixView<T>& A, const VectorView<T>& xD, 
      size_t* P, T& det, size_t j1=0)
  {
    // Bunch-Kauffman algorithm for LDL Decomposition of a symmetric matrix.
    //
    // We want to decompose the matrix (input as A) into P * L * D * Lt * Pt,
    // where P is a permutation, L is a lower triangle matrix with 1's for 
    // the diagonal, and D is a pseudo-diagonal matrix, meaning that 
    // there may be some 2x2 blocks along the diagonal.
    //
    // The reason LU decomposition is hard for symmetric matrices lies in
    // the pivoting.  If we pivot an off-diagonal element onto the diagonal,
    // then it destroys the symmetric structure of the matrix.
    // We can symmetrically pivot to move diagonal elements up or down the 
    // diagonal, but if a diagonal element is 0 (or very small compared to
    // an off-diagonal element) dividing by it causes problems.
    //
    // Bunch and Parlett showed that we can pivot off-diagonal elements 
    // to the first off diagonal and leave the 2x2 block unreduced in the 
    // D matrix.
    // When we are done, dividing by the 2x2 blocks is trivial to do
    // with correct pivoting.
    // 
    // (Following the explanation of Golub and van Loan, section 4.4.4:)
    // Consider the first step in the decomposition where we find
    // P1, E, C, B such that
    //
    // P1 A P1t = [ E Ct ]
    //            [ C B  ]
    // where E is sxs, C is (n-s)xs, B is (n-s)x(n-s), and s is either 1 or 2.
    // 
    // If A is non-zero, then E can be chosen to be non-singular, so:
    //
    // [ E Ct ] = [   I   0 ] [ E     0     ] [ I E^-1Ct ]
    // [ C B  ]   [ CE^-1 I ] [ 0 B-CE^-1Ct ] [ 0   I    ]
    //
    // CE^-1 is the first portion of the L matrix will will be making.
    // E is the first part of D.
    // B-CE^-1Ct = A~ is the submatrix on which to continue the process.
    //
    // Bunch and Parlett showed how to find the appropriate P,E to 
    // minimize the growth of A~.
    // Define mu0 = max_i,j |A(i,j)|
    //        mu1 = max_i |A(i,i)|
    //        
    // When the largest element is on the diagonal (mu0 = mu1), then 
    // it seems obvious that we would want to use it for E.
    // When the largest element is off-diagonal and it is much larger than
    // any on-diagonal element mu0 >> mu1, we want E to be 2x2,
    // with mu0 being the off-diagonal of E.
    //
    // When mu0 is only slightly larger than mu1, it is less clear which
    // strategy is better.
    //
    // The Bunch-Parlett strategy is to take:
    // s=1, E=mu1     when mu1 > alpha * mu0.
    // s=2, E_10=mu0  when mu1 < alpha * mu0.
    // where alpha is some parameter in [0,1].
    //
    // To determine what is a good value for alpha, find the bound on
    // the values of A~ in each case:
    // s=1: |a~_ij| < max(|B-CE^-1Ct|) < max(B) + max(CCt/mu1) 
    //              < mu0 + max(CCt)/(alpha*mu0) 
    //              < mu0 + mu0^2/(alpha*mu0)
    //              = mu0 ( 1 + 1/alpha )
    // s=2: |a~_ij| < max(|B-CE^-1Ct|) < max(B) + max(CE^-1Ct)
    //              < mu0 + max(|[ mu0 mu0 ] [ x  mu0 ]^-1 [ mu0 ]|)
    //                                       [ mu0 y  ]    [ mu0 ]
    //              [ max over (x,y) with -alpha*mu0 < x,y < alpha*mu0 ]
    //              = mu0 + max(|[ mu0 mu0 ] [ y  -mu0 ] [ mu0 ] / (xy-mu0^2)|)
    //                                       [ -mu0  x ] [ mu0 ] 
    //              = mu0 + max( |x+y-2mu0| mu0^2 / |xy-mu0^2| )
    //              [ Numerator is largest when x+y = -2alpha*mu0.        ]
    //              [ Denominator is smallest when xy = alpha^2 mu0^2.    ]
    //              [ These are consistent with each other, so that's it. ]
    //              = mu0 + 2(1+alpha)mu0/(1-alpha^2)
    //              = mu0 ( 1-alpha^2+2+2*alpha ) / (1-alpha^2)
    //              = mu0 ( 3 + 2alpha - alpha^2) / (1-alpha^2)
    //              = mu0 ( 3 - alpha ) / ( 1 - alpha )
    //
    // The optimal alpha is where the growth from 2 s=1 steps equals the 
    // growth from a single s=2 step:
    //
    // (1+1/a)^2 = (3-a)/(1-a)
    // 1+a-a^2-a^3 = 3a^2 - a^3
    // 1+a-4a^2 = 0
    // alpha = (1+sqrt(17))/8
    //
    // The only trouble with this algorithm is that finding mu0 requires
    // O(n^2) comparisons, which needs to be done O(n) times, so these become 
    // a significant fraction of the computation time.
    //
    // Bunch and Kauffman modified the algorithm slightly to require only
    // O(n) comparisons for each step.  The values to calculate are:
    // a00 = |A(0,0)|
    // ap0 = |A(p,0)| = max_i |A(i,0)| 
    // apq = |A(p,q)| = max_j!=p |A(p,j)|
    // app = |A(p,p)|
    //
    // Then their tests are:
    //
    // if a00 > alpha * ap0
    //   s = 1, E = a00  (No need to calculate arp.)
    // else if a00*apq > alpha * ap0^2
    //   s = 1, E = a00
    // else if app > alpha * apq
    //   s = 1, E = app
    // else 
    //   s = 2, E = (a00,app,ap0)
    //   [Note: Golub and van Loan wrongly say to put apq in E, 
    //          rather than ap0 here.]
    // 
    TMVAssert(A.size()>0);
    TMVAssert(xD.size()+1 == A.size());
    TMVAssert(IsReal(T()) || herm == A.isherm());
    TMVAssert(herm || IsComplex(T()));
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);
    const size_t N = A.size();
    static const RealType(T) alpha = (SQRT(RealType(T)(17))+1)/8;

#ifdef XDEBUG
    Matrix<T> L(N,N,T(0));
    LowerTriMatrixViewOf(L) = A.LowerTri(UnitDiag);
    L.SubMatrix(j1,N,j1,N).SetToIdentity();
    Matrix<T> DD(N,N,T(0));
    DD.diag() = A.diag();
    DD.diag(-1) = xD;
    DD.diag(1) = herm ? xD.Conjugate() : xD;
    DD.SubMatrix(j1,N,j1,N) = A.SubSymMatrix(j1,N);
    Matrix<T> A0 = L*DD*(herm ? L.Adjoint() : L.Transpose());
    for(size_t i=j1;i<A.size();i++) P[i]=i;
    A0.ReversePermuteRows(P);
    A0.ReversePermuteCols(P);
#endif

    VectorView<T> D = A.diag();

    for (size_t j=j1; j<N;) // ++j or j+=2 done below
    {
      bool seq1 = true;
      if (j == N-1) {
	// No permutation
	P[j] = j;
      } else {
	RealType(T) ajj = herm ? abs(REAL(D(j))) : abs(D(j));
	size_t p; // p is relative to j index, not absolute.
	RealType(T) apj = A.col(j,j,N).MaxAbsElement(&p);

	if (p == 0 || ajj >= alpha * apj) {
	  // No permutation
	  P[j] = j;
	} else {
	  p+=j;
	  RealType(T) app = herm ? abs(REAL(D(p))) : abs(D(p));
	  RealType(T) apq = A.row(p,j,p).MaxAbsElement();
	  if (p+1 < N) {
	    RealType(T) apq2 = A.col(p,p+1,N).MaxAbsElement();
	    apq = max(apq,apq2);
	  }
	  if (ajj*apq >= alpha * apj * apj) {
	    // No permutation
	    P[j] = j;
	  } else if (app >= alpha * apq) {
	    // Permute p diagonal into j spot
	    TMVAssert(p<A.size());
	    A.SwapRowsCols(j,p);
	    P[j] = p;
#ifdef XDEBUG
	    L.SubMatrix(0,N,0,j).SwapRows(j,p);
#endif
	  } else {
	    // Permute pj element into j+1,j spot
	    // This also permutes pp element into j+1,j+1
	    seq1 = false;
	    P[j] = j;
	    P[j+1] = p;
	    if (p != j+1) {
	      TMVAssert(p<A.size());
	      A.SwapRowsCols(j+1,p);
#ifdef XDEBUG
	      L.SubMatrix(0,N,0,j).SwapRows(j+1,p);
#endif
	    }
	  }
	}
      }

      // Now the LU solving:
      if (seq1) {
	if (herm) {
	  RealType(T) dj = REAL(D(j));
	  det *= dj;
	  if (dj != RealType(T)(0)) {
	    A.col(j,j+1,N) /= dj;
	    A.SubSymMatrix(j+1,N) -= 
	      dj*A.col(j,j+1,N)^A.col(j,j+1,N).Conjugate();
	  }
	} else {
	  T dj = D(j);
	  det *= dj;
	  if (dj != T(0)) {
	    A.col(j,j+1,N) /= dj;
	    A.SubSymMatrix(j+1,N) -= dj*A.col(j,j+1,N)^A.col(j,j+1,N);
	  }
	}
#ifdef XDEBUG
	DD.col(j,j+1,N).Zero();
	DD.row(j,j+1,N).Zero();
	DD(j,j) = D(j);
	L.col(j,j+1,N) = A.col(j,j+1,N);
	DD.SubMatrix(j+1,N,j+1,N) = A.SubSymMatrix(j+1,N);
        Matrix<T> A2 = L*DD*(herm?L.Adjoint():L.Transpose());
	A2.ReversePermuteRows(P);
	A2.ReversePermuteCols(P);
	if (Norm(A0-A2) > 1.e-5*Norm(A0)) {
	  if (herm)
	    cerr<<"Herm: s==1\n";
	  else
	    cerr<<"Sym: s==1\n";
	  cerr<<"j = "<<j<<endl;
	  cerr<<"A = "<<A<<endl;
	  cerr<<"D = "<<D<<endl;
	  cerr<<"xD = "<<xD<<endl;
	  cerr<<"L = "<<L<<endl;
	  cerr<<"DD = "<<DD<<endl;
	  cerr<<"A2 = "<<A2<<endl;
	  cerr<<"A0 = "<<Type(A)<<"  "<<A0<<endl;
	  cerr<<"Norm(A0-A2) = "<<Norm(A0-A2)<<endl;
	  abort();
	}
#endif
	++j;
      } else {
	// Invert E:  E^-1 = [ x z* ]
	//                   [ z y  ]

	T x = D(j);
	T y = D(j+1);
	T z = xD(j) = A(j+1,j);
	A(j+1,j)=T(0);
	T d;
	SymInvert_2x2<herm>(x,y,z,&d);
	if (herm) det *= REAL(d);
	else det *= d;

	if (A.iscm()) {
	  // Call C the current kx2 matrix which is equal to LE
	  // Need to right-multiply by E^-1 to get L
	  // A(j+2:N,j:j+2) = CE^-1
	  //
	  // Also need to update the rest of the matrix by -= L * Ct
	  // A(j+2:N,j+2:N) -= CE^-1Ct 
	  //                -= C (CE^-1)t (since E^-1 is Hermitian)
	  // A(p,q) -= C(p,0)*(L(q,0)*) + C(p,1)*(L(q,1)*)
	  //   p >= q, so loop from top down
	  const int sj = A.stepj();
	  T* Cq0 = A.ptr() + (j+2) + j*sj;
	  T* Cq1 = Cq0 + sj;
	  for(size_t q = j+2; q<N; ++q,++Cq0,++Cq1) {
	    T Lq0 = *Cq0 * x + *Cq1 * z;
	    T Lq1 = *Cq0 * (herm?CONJ(z):z) + *Cq1 * y;
	    A.col(q,q,N) -= A.col(j,q,N) * (herm?CONJ(Lq0):Lq0);
	    A.col(q,q,N) -= A.col(j+1,q,N) * (herm?CONJ(Lq1):Lq1);
	    *Cq0 = Lq0;
	    *Cq1 = Lq1;
	  }
	} else {
	  // A(j+2:N,j:j+2) = CE^-1
	  // A(j+2:N,j+2:N) -= CE^-1Ct
	  // A(p,q) -= L(p,0)*(C(q,0)*) + L(p,1)*(C(q,1)*)
	  //   p >= q, so loop from bottom up.
	  const int si = A.stepi();
	  T* Cp0 = A.ptr() + (N-1)*si + j;
	  T* Cp1 = Cp0 + 1;
	  for(size_t p = N-1; p>=j+2; --p,Cp0-=si,Cp1-=si) {
	    T Lp0 = *Cp0 * x + *Cp1 * z;
	    T Lp1 = *Cp0 * (herm?CONJ(z):z) + *Cp1 * y;
	    A.row(p,j+2,p+1) -= Lp0 * 
	      (herm ? A.col(j,j+2,p+1).Conjugate() : A.col(j,j+2,p+1));
	    A.row(p,j+2,p+1) -= Lp1 *
	      (herm ? A.col(j+1,j+2,p+1).Conjugate() : A.col(j+1,j+2,p+1));
	    *Cp0 = Lp0;
	    *Cp1 = Lp1;
	  }
	}

#ifdef XDEBUG
	DD.col(j,j+2,N).Zero();
	DD.col(j+1,j+2,N).Zero();
	DD.row(j,j+2,N).Zero();
	DD.row(j+1,j+2,N).Zero();
	DD(j,j) = D(j);
	DD(j+1,j+1) = D(j+1);
	DD(j+1,j) = xD(j);
	DD(j,j+1) = herm ? CONJ(xD(j)) : xD(j);
	L.col(j,j+2,N) = A.col(j,j+2,N);
	L.col(j+1,j+2,N) = A.col(j+1,j+2,N);
	DD.SubMatrix(j+2,N,j+2,N) = A.SubSymMatrix(j+2,N);
        Matrix<T> A2 = L*DD*(herm?L.Adjoint():L.Transpose());
	A2.ReversePermuteRows(P);
	A2.ReversePermuteCols(P);
	if (Norm(A0-A2) > 1.e-5*Norm(A0)) {
	  if (herm)
	    cerr<<"Herm: s==2\n";
	  else
	    cerr<<"Sym: s==2\n";
	  cerr<<"j = "<<j<<endl;
	  cerr<<"A = "<<A<<endl;
	  cerr<<"D = "<<D<<endl;
	  cerr<<"xD = "<<xD<<endl;
	  cerr<<"L = "<<L<<endl;
	  cerr<<"DD = "<<DD<<endl;
	  cerr<<"A2 = "<<A2<<endl;
	  cerr<<"A0 = "<<Type(A)<<"  "<<A0<<endl;
	  cerr<<"Norm(A0-A2) = "<<Norm(A0-A2)<<endl;
	  abort();
	}
#endif
	j+=2;
      }
    }
#ifdef XDEBUG
    DD = DiagMatrixViewOf(D);
    DD.diag(-1) = xD;
    DD.diag(1) = herm?xD.Conjugate():xD;
    L = A.LowerTri(UnitDiag);
    Matrix<T> A2 = L*DD*(herm?L.Adjoint():L.Transpose());
    A2.ReversePermuteRows(P);
    A2.ReversePermuteCols(P);

    if (Norm(A2-A0) > 1.e-5*SQR(Norm(L))*Norm(DD)) {
      cerr<<"NonBlockSymLUDecompose\n";
      cerr<<"A0 = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"D = "<<D<<endl;
      cerr<<"xD = "<<xD<<endl;
      cerr<<"L = "<<L<<endl;
      cerr<<"DD = "<<DD<<endl;
      cerr<<"L*D*Lt = "<<L*DD*L.Adjoint()<<endl;
      cerr<<"A2 = "<<A2<<endl;
      cerr<<"A2-A0 = "<<Matrix<T>(A2-A0).Clip(1.e-5)<<endl;
      cerr<<"Norm(A0-A2) = "<<Norm(A0-A2)<<endl;
      cerr<<"nm = "<<Norm(A2-A0)/(SQR(Norm(L))*Norm(DD))<<endl;
      abort();
    }
#endif
  }

  template <bool herm, class T> inline void BlockSymLU_Decompose(
      const SymMatrixView<T>& A, const VectorView<T>& xD, 
      size_t* P, T& det)
  {
    // Blocked version with level 3 calls

    TMVAssert(A.size()>0);
    TMVAssert(xD.size()+1 == A.size());
    TMVAssert(IsReal(T()) || herm == A.isherm());
    TMVAssert(herm || IsComplex(T()));
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);

#ifdef XDEBUG
    Matrix<T> A0 = A;
#endif

    static const RealType(T) alpha = (SQRT(RealType(T)(17))+1)/8;
    const size_t N = A.size();

    VectorView<T> D = A.diag();
    // Only make LD if we're going to use it:
    Matrix<T,ColMajor> LD(
	N < SYM_LU_BLOCKSIZE ? 0 : N,
	N < SYM_LU_BLOCKSIZE ? 0 : SYM_LU_BLOCKSIZE+1);

    for (size_t j1=0; j1<N; ) {
      size_t j2 = min(j1+SYM_LU_BLOCKSIZE,N);

      if (j2 < N) { // On last loop we skip some steps.  See below.
	size_t j;
	for (j=j1; j<j2;) { // ++j or j+=2 done below
	  bool seq1 = true;
	  size_t jj = j-j1; // Col index in LD

	  LD.col(jj,j,N) = A.col(j,j,N);
	  LD.col(jj,j,N) -= A.SubMatrix(j,N,j1,j) *
	    (herm ? LD.row(j,0,jj).Conjugate() : LD.row(j,0,jj));

	  if (j == N-1) {
	    // No permutation
	    P[j] = j;
	  } else {
	    RealType(T) ajj = abs(LD(j,jj));
	    size_t p; // p is relative to j+1 index, not absolute.
	    RealType(T) apj = LD.col(jj,j+1,N).MaxAbsElement(&p);

	    if (ajj >= alpha * apj) {
	      // No permutation
	      P[j] = j;
	    } else {
	      p+=j+1;

	      LD.col(jj+1,j,p) = A.col(p,j,p);
	      LD.col(jj+1,p,N) = A.col(p,p,N);
	      LD.col(jj+1,j,N) -= A.SubMatrix(j,N,j1,j) * 
		(herm ? LD.row(p,0,jj).Conjugate() : LD.row(p,0,jj));

	      RealType(T) app = abs(LD(p,jj+1));
	      RealType(T) apq = LD.col(jj+1,j,p).MaxAbsElement();
	      if (p+1 < N) {
		RealType(T) apq2 = LD.col(jj+1,p+1,N).MaxAbsElement();
		apq = max(apq,apq2);
	      }
	      if (ajj*apq >= alpha * apj * apj) {
		// No permutation
		P[j] = j;
	      } else if (app >= alpha * apq) {
		// Permute p diagonal into j spot
		P[j] = p;
		// Move updated pivot column (in LD) to j column
		LD.col(jj,j,N) = LD.col(jj+1,j,N);
		TMVAssert(p<A.size());
		// Do important parts of the A.SwapRowsCols(j,p) call,
		// We will overwrite A.col(j), so don't bother swapping into it.
		D(p) = D(j);
		A.row(p,j+1,p) = 
		  herm ? A.col(j,j+1,p).Conjugate() : A.col(j,j+1,p);
		A.col(p,p+1,N) = A.col(j,p+1,N);
		Swap(A.row(j,0,j),A.row(p,0,j));
		// Note: this swap goes all the way back to 0, not j1
		// Also need to do row swaps in LD
		Swap(LD.row(j,0,jj+1),LD.row(p,0,jj+1));
	      } else {
		// Permute pj element into j+1,j spot
		// This also permutes pp element into j+1,j+1
		seq1 = false;
		P[j] = j;
		P[j+1] = p;
		if (p != j+1) {
		  TMVAssert(p<A.size());
		  // Do important parts of the A.SwapRowsCols(j+1,p) call,
		  // We will overwrite A.cols(j,j+1), so don't bother 
		  // swapping into them.
		  D(p) = D(j+1);
		  A.row(p,j+2,p) = 
		    (herm ? A.col(j+1,j+2,p).Conjugate() : A.col(j+1,j+2,p));
		  A.col(p,p+1,N) = A.col(j+1,p+1,N);
		  Swap(A.row(j+1,0,j),A.row(p,0,j));
		  // Also need to do row swaps in LD
		  Swap(LD.row(j+1,0,jj+2),LD.row(p,0,jj+2));
		}
	      }
	    }
	  }

	  // Now the LU solving:
	  if (seq1) {
	    // LD.col(j) now holds L.col(j)*D(j)
	    A.col(j,j,N) = LD.col(jj,j,N);
	    if (herm) {
	      RealType(T) dj = REAL(D(j));
	      det *= dj;
	      if (dj != RealType(T)(0)) A.col(j,j+1,N) /= dj;
	    } else {
	      T dj = D(j);
	      det *= dj;
	      if (dj != T(0)) A.col(j,j+1,N) /= dj;
	    }
	    ++j;
	  } else {
	    // LD.cols(j,j+1) now hold L.cols(j,j+1) * E 
	    // Invert E:  E^-1 = [ x z* ]
	    //                   [ z y  ]

	    if (herm) {
	      D(j) = REAL(LD(j,jj));
	      D(j+1) = REAL(LD(j+1,jj+1));
	    } else {
	      D(j) = LD(j,jj);
	      D(j+1) = LD(j+1,jj+1);
	    }
	    T x = D(j);
	    T y = D(j+1);
	    T z = xD(j) = LD(j+1,jj);
	    A(j+1,j)=T(0);

	    T d;
	    SymInvert_2x2<herm>(x,y,z,&d);
	    if (herm) det *= REAL(d);
	    else det *= d;

	    // A(j+2:N,j:j+2) = LD E^-1
	    const int ldsj = LD.stepj();
	    T* LDq0 = LD.ptr() + (j+2) + jj*ldsj;
	    T* LDq1 = LDq0 + ldsj;
	    if (A.iscm()) {
	      const int sj = A.stepj();
	      T* Aq0 = A.ptr() + (j+2) + j*sj;
	      T* Aq1 = Aq0 + sj;
	      for(size_t q = j+2; q<N; ++q,++Aq0,++Aq1,++LDq0,++LDq1) {
		*Aq0 = *LDq0 * x + *LDq1 * z;
		*Aq1 = *LDq0 * (herm?CONJ(z):z) + *LDq1 * y;
	      }
	    } else {
	      const int si = A.stepi();
	      T* Aq0 = A.ptr() + (j+2)*si + j;
	      T* Aq1 = Aq0 + 1;
	      for(size_t q = j+2; q<N; ++q,Aq0+=si,Aq1+=si,++LDq0,++LDq1) {
		*Aq0 = *LDq0 * x + *LDq1 * z;
		*Aq1 = *LDq0 * (herm?CONJ(z):z) + *LDq1 * y;
	      }
	    }
	    j+=2;
	  }
	}
	j2 = j; // in case last one was a j+=2

	// Now update the rest of A:
	// A(j2:N,j2:N) -= L(j2:N,j1:j2)*D(j1:j2,j1:j2)*L(j2:N,j1:j2)t
	// A(j2:N,j2:N) -= L(j2:N,j1:j2)*LD(j2:N,j1:j2)t
	SymMatrixView<T> A22 = A.SubSymMatrix(j2,N);
	MatrixView<T> L21 = A.SubMatrix(j2,N,j1,j2);
	MatrixView<T> LD21 = LD.SubMatrix(j2,N,0,j2-j1);

	// A22 -= L21 * LD21t
	SymMultMM(T(-1),L21,herm?LD21.Adjoint():LD21.Transpose(),1,A22);

      } else NonBlockSymLU_Decompose<herm>(A,xD,P,det,j1);

#ifdef XDEBUG
      Matrix<T> L = A.LowerTri();
      L.Cols(j2,N).Zero();
      L.diag().SetAllTo(T(1));

      Matrix<T> DD = DiagMatrixViewOf(A.diag());
      DD.diag(-1) = xD;
      DD.diag(1) = herm?xD.Conjugate():xD;
      DD.SubMatrix(j2,N,j2,N) = A.SubSymMatrix(j2,N);

      Matrix<T> A2 = L * DD * (herm?L.Adjoint():L.Transpose());
      for(size_t j=j2;j<N;j++) P[j] = j;
      A2.ReversePermuteRows(P);
      A2.ReversePermuteCols(P);

      if (Norm(A2-A0) > 1.e-5*SQR(Norm(L))*Norm(DD)) {
	cerr<<"BlockSymLUDecompose\n";
	cerr<<"j1 = "<<j1<<", j2 = "<<j2<<", nb = "<<SYM_LU_BLOCKSIZE<<endl;
	cerr<<"A0 = "<<Type(A)<<A0<<endl;

	Vector<T> xD1(xD.size(),T(0));
	size_t P1[A.size()];
	T d1(1);;
	Matrix<T,ColMajor> m1(A.size(),A.size());
	SymMatrixView<T> A1 = herm ?
	  HermMatrixViewOf(m1,Lower) : SymMatrixViewOf(m1,Lower);
	A1.LowerTri() = LowerTriMatrixViewOf(A0);
	cerr<<"A1 = "<<Type(A1)<<A1<<endl;
	NonBlockSymLU_Decompose<herm>(A1.View(),xD1.View(),P1,d1);
	cerr<<"L1 = "<<A1.LowerTri(UnitDiag)<<endl;
	cerr<<"D1 = "<<A1.diag()<<endl;
	cerr<<"xD1 = "<<xD1<<endl;
	cerr<<"P1 = ";
	for(size_t i=0;i<A.size();i++) cerr<<P1[i]<<" ";
	cerr<<endl;

	cerr<<"L = "<<L<<endl;
	cerr<<"D = "<<A.diag().Real()<<endl;
	cerr<<"xD = "<<xD<<endl;
	cerr<<"P = ";
	for(size_t i=0;i<A.size();i++) cerr<<P[i]<<" ";
	cerr<<endl;

	cerr<<"L*D*Lt = "<<L*DD*(herm?L.Adjoint():L.Transpose())<<endl;
	cerr<<"A2 = "<<A2<<endl;
	cerr<<"A2-A0 = "<<A2-A0<<endl;
	cerr<<"nm = "<<Norm(A2-A0)/(SQR(Norm(L))*Norm(DD))<<endl;

	abort();
      }
#endif
      j1 = j2;
    }
  }

  template <bool herm, class T> inline void NonLapSymLU_Decompose(
      const SymMatrixView<T>& A, const VectorView<T>& xD, 
      size_t* P, T& det)
  {
    TMVAssert(A.size()>0);
    TMVAssert(xD.size()+1 == A.size());
    TMVAssert(herm == A.isherm());
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);

    if (A.size() > SYM_LU_BLOCKSIZE) 
      BlockSymLU_Decompose<herm>(A,xD,P,det);
    else 
      NonBlockSymLU_Decompose<herm>(A,xD,P,det);
  }

#ifdef LAP
  template <class T> inline void LapSymLU_Decompose(
      const SymMatrixView<T>& A, const VectorView<T>& xD,
      size_t* P, T& det)
  { NonLapSymLU_Decompose<false>(A,xD,P,det); }
  template <> inline void LapSymLU_Decompose(
      const SymMatrixView<double>& A, const VectorView<double>& xD,
      size_t* P, double& det)
  {
    TMVAssert(A.size()>0);
    TMVAssert(xD.size()+1 == A.size());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);
    TMVAssert(A.iscm());

    int n = A.size();
    int lda = A.stepj();
    auto_array<int> lap_p(new int[n]);
#ifndef LAPNOWORK
    int lwork = n*LAP_BLOCKSIZE;
    double* work = LAP_DWork(lwork);
#endif
    LAPNAME(dsytrf) (LAPCM LAPCH_LO,LAPV(n),LAPP(A.ptr()),LAPV(lda),
	LAPP(lap_p.get()) LAPWK(work) LAPVWK(lwork) LAPINFO LAP1);
    LAP_Results("dsytrf");
    int* pi = lap_p.get();
    for(size_t i=0;i<A.size();++i,++pi) {
      int ip = *pi LAPMINUS1;
      if (ip != int(i)) {
	if (ip >= 0) {
	  Swap(A.row(i,0,i),A.row(ip,0,i));
	  P[i] = ip;
	  det *= A(i,i);
	} else {
	  TMVAssert(i+1 < A.size());
	  TMVAssert(*(pi+1) == *pi);
	  ip = -*pi LAPMINUS1;
	  P[i] = i;
	  P[i+1] = ip;
	  if (int(i+1) != ip) {
	    Swap(A.row(i+1,0,i),A.row(ip,0,i));
	  }
	  double& xDi = A(i+1,i);
	  xD(i) = xDi;
	  det *= A(i,i)*A(i+1,i+1)-xDi*xDi;
	  xDi = double(0);
	  ++i; ++pi; // extra ++ for 2x2 case
	}
      } else {
	det *= A(i,i);
	P[i] = i;
      }
    }
  }
  template <> inline void LapSymLU_Decompose(
      const SymMatrixView<complex<double> >& A,
      const VectorView<complex<double> >& xD,
      size_t* P, complex<double>& det)
  {
    TMVAssert(A.size()>0);
    TMVAssert(xD.size()+1 == A.size());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);
    TMVAssert(A.iscm());
    int n = A.size();
    int lda = A.stepj();
    auto_array<int> lap_p(new int[n]);
#ifndef LAPNOWORK
    int lwork = n*LAP_BLOCKSIZE;
    complex<double>* work = LAP_ZWork(lwork);
#endif
    if (A.isherm()) {
      LAPNAME(zhetrf) (LAPCM LAPCH_LO,LAPV(n),LAPP(A.ptr()),LAPV(lda),
	  LAPP(lap_p.get()) LAPWK(work) LAPVWK(lwork) LAPINFO LAP1);
      LAP_Results("zhetrf");
    } else {
      LAPNAME(zsytrf) (LAPCM LAPCH_LO,LAPV(n),LAPP(A.ptr()),LAPV(lda),
	  LAPP(lap_p.get()) LAPWK(work) LAPVWK(lwork) LAPINFO LAP1);
      LAP_Results("zsytrf");
    }
    int* pi = lap_p.get();
    for(size_t i=0;i<A.size();++i,++pi) {
      int ip = *pi LAPMINUS1;
      if (ip != int(i)) {
	if (ip >= 0) {
	  Swap(A.row(i,0,i),A.row(ip,0,i));
	  P[i] = ip;
	  if (A.isherm()) det *= REAL(A(i,i));
	  else det *= A(i,i);
	} else {
	  TMVAssert(i+1 < A.size());
	  TMVAssert(*(pi+1) == *pi);
	  ip = -*pi LAPMINUS1;
	  P[i] = i;
	  P[i+1] = ip;
	  if (int(i+1) != ip) {
	    Swap(A.row(i+1,0,i),A.row(ip,0,i));
	  }
	  complex<double>& xDi = A(i+1,i).GetRef();
	  xD(i) = xDi;
	  if (A.isherm()) det *= REAL(A(i,i))*REAL(A(i+1,i+1))-NORM(xDi);
	  else det *= A(i,i)*A(i+1,i+1)-xDi*xDi;
	  xDi = double(0);
	  ++i; ++pi; // extra ++ for 2x2 case
	}
      } else {
	if (A.isherm()) det *= REAL(A(i,i));
	else det *= A(i,i);
	P[i] = i;
      }
    }
  }
#ifndef NOFLOAT
  template <> inline void LapSymLU_Decompose(
      const SymMatrixView<float>& A, const VectorView<float>& xD,
      size_t* P, float& det)
  {
    TMVAssert(A.size()>0);
    TMVAssert(xD.size()+1 == A.size());
    TMVAssert(A.isherm());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);
    TMVAssert(A.iscm());
    TMVAssert(A.isherm());
    int n = A.size();
    int lda = A.stepj();
    auto_array<int> lap_p(new int[n]);
#ifndef LAPNOWORK
    int lwork = n*LAP_BLOCKSIZE;
    float* work = LAP_SWork(lwork);
#endif
    LAPNAME(ssytrf) (LAPCM LAPCH_LO,LAPV(n),LAPP(A.ptr()),LAPV(lda),
	LAPP(lap_p.get()) LAPWK(work) LAPVWK(lwork) LAPINFO LAP1);
    LAP_Results("ssytrf");

    int* pi = lap_p.get();
    for(size_t i=0;i<A.size();++i,++pi) {
      int ip = *pi LAPMINUS1;
      if (ip != int(i)) {
	if (ip >= 0) {
	  Swap(A.row(i,0,i),A.row(ip,0,i));
	  P[i] = ip;
	  det *= A(i,i);
	} else {
	  TMVAssert(i+1 < A.size());
	  TMVAssert(*(pi+1) == *pi);
	  ip = -*pi LAPMINUS1;
	  P[i] = i;
	  P[i+1] = ip;
	  if (int(i+1) != ip) {
	    Swap(A.row(i+1,0,i),A.row(ip,0,i));
	  }
	  float& xDi = A(i+1,i);
	  xD(i) = xDi;
	  det *= A(i,i)*A(i+1,i+1)-xDi*xDi;
	  xDi = float(0);
	  ++i; ++pi; // extra ++ for 2x2 case
	}
      } else {
	det *= A(i,i);
	P[i] = i;
      }
    }
  }
  template <> inline void LapSymLU_Decompose(
      const SymMatrixView<complex<float> >& A,
      const VectorView<complex<float> >& xD,
      size_t* P, complex<float>& det)
  {
    TMVAssert(A.size()>0);
    TMVAssert(xD.size()+1 == A.size());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);
    TMVAssert(A.iscm());
    int n = A.size();
    int lda = A.stepj();
    auto_array<int> lap_p(new int[n]);
#ifndef LAPNOWORK
    int lwork = n*LAP_BLOCKSIZE;
    complex<float>* work = LAP_CWork(lwork);
#endif
    if (A.isherm()) {
      LAPNAME(chetrf) (LAPCM LAPCH_LO,LAPV(n),LAPP(A.ptr()),LAPV(lda),
	  LAPP(lap_p.get()) LAPWK(work) LAPVWK(lwork) LAPINFO LAP1);
      LAP_Results("chetrf");
    } else {
      LAPNAME(csytrf) (LAPCM LAPCH_LO,LAPV(n),LAPP(A.ptr()),LAPV(lda),
	  LAPP(lap_p.get()) LAPWK(work) LAPVWK(lwork) LAPINFO LAP1);
      LAP_Results("csytrf");
    }
    int* pi = lap_p.get();
    for(size_t i=0;i<A.size();++i,++pi) {
      int ip = *pi LAPMINUS1;
      if (ip != int(i)) {
	if (ip >= 0) {
	  Swap(A.row(i,0,i),A.row(ip,0,i));
	  P[i] = ip;
	  if (A.isherm()) det *= REAL(A(i,i));
	  else det *= A(i,i);
	} else {
	  TMVAssert(i+1 < A.size());
	  TMVAssert(*(pi+1) == *pi);
	  ip = -*pi LAPMINUS1;
	  P[i] = i;
	  P[i+1] = ip;
	  if (int(i+1) != ip) {
	    Swap(A.row(i+1,0,i),A.row(ip,0,i));
	  }
	  complex<float>& xDi = A(i+1,i).GetRef();
	  xD(i) = xDi;
	  if (A.isherm()) det *= REAL(A(i,i))*REAL(A(i+1,i+1))-NORM(xDi);
	  else det *= A(i,i)*A(i+1,i+1)-xDi*xDi;
	  xDi = float(0);
	  ++i; ++pi; // extra ++ for 2x2 case
	}
      } else {
	if (A.isherm()) det *= REAL(A(i,i));
	else det *= A(i,i);
	P[i] = i;
      }
    }
  }
#endif // NOFLOAT
#endif // LAP

  template <class T> void SymLU_Decompose(
      const SymMatrixView<T>& A, const VectorView<T>& xD, 
      size_t* P, T& det)
  {
    TMVAssert(A.size() > 0);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(xD.size()+1 == A.size());

#ifdef XDEBUG
    Matrix<T> A0 = A;
#endif

    if (A.isconj()) SymLU_Decompose(A.Conjugate(),xD,P,det);
    else if (A.uplo() == Upper) {
      if (A.isherm()) SymLU_Decompose(A.Adjoint(),xD,P,det);
      else SymLU_Decompose(A.Transpose(),xD,P,det);
    }
    else {
      xD.Zero();
#ifdef LAP
      if (!A.iscm()) {
	Matrix<T,ColMajor> temp(A.size(),A.size());
	SymMatrixView<T> A2 = A.isherm() ?
	  HermMatrixViewOf(temp,Lower) :
	  SymMatrixViewOf(temp,Lower);
	A2 = A;
	LapSymLU_Decompose(A2,xD,P,det);
	A = A2;
      } else {
	LapSymLU_Decompose(A,xD,P,det);
      }
#else
      if (A.isherm()) NonLapSymLU_Decompose<true>(A,xD,P,det);
      else NonLapSymLU_Decompose<false>(A,xD,P,det);
#endif
    }

#ifdef XDEBUG
    LowerTriMatrix<T,UnitDiag> L = A.LowerTri(UnitDiag);
    Matrix<T> DD(A.size(),A.size(),T(0));
    DD.diag() = A.diag();
    DD.diag(-1) = xD;
    DD.diag(1) = A.isherm() ? xD.Conjugate() : xD;
    Matrix<T> A2 = L*DD*(A.isherm() ? L.Adjoint() : L.Transpose());
    A2.ReversePermuteRows(P);
    A2.ReversePermuteCols(P);
    if (Norm(A2-A0) > 0.001*Norm(A0)*NormSq(L)*Norm(DD)) {
      cerr<<"SymLU_Decompose\n";
      cerr<<"A0 = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"A -> "<<A<<endl;
      cerr<<"L = "<<L<<endl;
      cerr<<"D = "<<DD<<endl;
      cerr<<"A2 = "<<A2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_SymLUDiv_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


