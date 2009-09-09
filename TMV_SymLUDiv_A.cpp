
#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Sym.h"
#include "TMV_Diag.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
  const size_t SYM_LU_BLOCKSIZE = TMV_BLOCKSIZE;
  const size_t SYM_LU_BLOCKSIZE1 = TMV_BLOCKSIZE;
  const size_t SYM_LU_BLOCKSIZE2 = TMV_BLOCKSIZE/2;
#else
  const size_t SYM_LU_BLOCKSIZE = 48;
  const size_t SYM_LU_BLOCKSIZE1 = 64;
  const size_t SYM_LU_BLOCKSIZE2 = 32;
#endif

  //
  // Decompose
  //


  template <class T> void SymInvert_2x2(T& a, T& b, T& c, T* dd=0)
  {
    // Invert matrix [ a  c ] -->  1/(ab-c^2) [ b -c ]
    //               [ c  b ]                 [ -c a ]
    T d = a*b-c*c;
    swap(a,b);
    a/=d;
    b/=d;
    c = -c/d;
    if (dd) *dd = d;
  }

  template <class T> void HermInvert_2x2(
      RealType(T)& a, RealType(T)& b, T& c, RealType(T)* dd=0)
  {
    // Invert matrix [ a  c* ] -->  1/(ab-|c|^2) [ b -c* ]
    //               [ c  b  ]                   [ -c a  ]
    RealType(T) d = a*b-NORM(c);
    swap(a,b);
    a/=d;
    b/=d;
    c = -c/d;
    if (dd) *dd = d;
  }

  // MJ: This is a bit slower for row major, even with blocking, so write the
  // row major version where A is decomposed into Lt D L rather than L D Lt.
  template <class T> void NonBlockHermLU_Decompose(
      const SymMatrixView<T>& A, const VectorView<T>& xD, 
      size_t* P, RealType(T)& det, size_t j1=0)
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
    TMVAssert(A.isherm());
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);
    const size_t N = A.size();
    static const RealType(T) alpha = (SQRT(RealType(T)(17))+1)/8;

    VectorView<RealType(T)> D = A.diag().Real();
#ifdef XDEBUG
    Matrix<T> L(N,N,T(0));
    LowerTriMatrixViewOf(L) = A.LowerTri().MakeUnitDiag();
    L.SubMatrix(j1,N,j1,N).SetToIdentity();
    Matrix<T> DD(N,N,T(0));
    DD.diag() = A.diag();
    DD.diag(-1) = xD;
    DD.diag(1) = xD.Conjugate();
    DD.SubMatrix(j1,N,j1,N) = A.SubSymMatrix(j1,N);
    Matrix<T> A0 = L*DD*L.Adjoint();
    for(size_t i=j1;i<A.size();i++) P[i]=i;
    A0.ReversePermuteRows(P);
    A0.ReversePermuteCols(P);
#endif

    for (size_t j=j1; j<N;) // ++j or j+=2 done below
    {
      bool seq1 = true;
      if (j == N-1) {
	// No permutation
	P[j] = j;
      } else {
	RealType(T) ajj = abs(D(j));
	size_t p; // p is relative to j index, not absolute.
	RealType(T) apj = A.col(j,j,N).MaxAbsElement(&p);

	if (p == 0 || ajj >= alpha * apj) {
	  // No permutation
	  P[j] = j;
	} else {
	  p+=j;
	  RealType(T) app = abs(D(p));
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
	RealType(T) dj = D(j);
	if (dj == T(0))
	  tmv_error("Zero pivot found in HermLU_Decompose");
	det *= dj;
	A.col(j,j+1,N) /= dj;
	A.SubSymMatrix(j+1,N) -= 
	    dj*A.col(j,j+1,N)^A.col(j,j+1,N).Conjugate();
#ifdef XDEBUG
	DD.col(j,j+1,N).Zero();
	DD.row(j,j+1,N).Zero();
	DD(j,j) = D(j);
	L.col(j,j+1,N) = A.col(j,j+1,N);
	DD.SubMatrix(j+1,N,j+1,N) = A.SubSymMatrix(j+1,N);
        Matrix<T> A2 = L*DD*L.Adjoint();
	A2.ReversePermuteRows(P);
	A2.ReversePermuteCols(P);
	if (Norm(A0-A2) > 1.e-5*max(RealType(T)(1),Norm(A0))) {
	  cerr<<"Herm: s==1\n";
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
	RealType(T) x = D(j);
	RealType(T) y = D(j+1);
	T z = xD(j) = A(j+1,j);

	A(j+1,j)=T(0);
	RealType(T) d;
	HermInvert_2x2(x,y,z,&d);
	det *= d;

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
	    T Lq1 = *Cq0 * CONJ(z) + *Cq1 * y;
	    A.col(q,q,N) -= A.col(j,q,N) * CONJ(Lq0);;
	    A.col(q,q,N) -= A.col(j+1,q,N) * CONJ(Lq1);;
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
	    T Lp1 = *Cp0 * CONJ(z) + *Cp1 * y;
	    A.row(p,j+2,p+1) -= Lp0 * A.col(j,j+2,p+1).Conjugate();
	    A.row(p,j+2,p+1) -= Lp1 * A.col(j+1,j+2,p+1).Conjugate();
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
	DD(j,j+1) = CONJ(xD(j));
	L.col(j,j+2,N) = A.col(j,j+2,N);
	L.col(j+1,j+2,N) = A.col(j+1,j+2,N);
	DD.SubMatrix(j+2,N,j+2,N) = A.SubSymMatrix(j+2,N);
        Matrix<T> A2 = L*DD*L.Adjoint();
	A2.ReversePermuteRows(P);
	A2.ReversePermuteCols(P);
	if (Norm(A0-A2) > 1.e-5*max(RealType(T)(1),Norm(A0))) {
	  cerr<<"Herm: s==2\n";
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
    Matrix<T> DD2 = DiagMatrixViewOf(D);
    DD2.diag(1) = xD.Conjugate();
    DD2.diag(-1) = xD;
    Matrix<T> L2 = A.LowerTri().MakeUnitDiag();
    Matrix<T> A2 = L2*DD2*L2.Adjoint();
    A2.ReversePermuteRows(P);
    A2.ReversePermuteCols(P);

    if (Norm(A2-A0) > 1.e-5*max(RealType(T)(1),SQR(Norm(L))*Norm(DD))) {
      cerr<<"NonBlockHermLUDecompose\n";
      cerr<<"A0 = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"D = "<<D<<endl;
      cerr<<"xD = "<<xD<<endl;
      cerr<<"L = "<<L<<endl;
      cerr<<"L2 = "<<L2<<endl;
      cerr<<"DD = "<<DD<<endl;
      cerr<<"DD2 = "<<DD2<<endl;
      cerr<<"L*D*Lt = "<<L*DD*L.Adjoint()<<endl;
      cerr<<"L2*D2*L2t = "<<L2*DD2*L2.Adjoint()<<endl;
      cerr<<"A2 = "<<A2<<endl;
      cerr<<"A2-A0 = "<<Matrix<T>(A2-A0).Clip(1.e-5)<<endl;
      cerr<<"Norm(A0-A2) = "<<Norm(A0-A2)<<endl;
      cerr<<"nm = "<<Norm(A2-A0)/(SQR(Norm(L))*Norm(DD))<<endl;
      abort();
    }
#endif
  }

  template <class T> void NonBlockSymLU_Decompose(
      const SymMatrixView<T>& A, const VectorView<T>& xD, 
      size_t* P, T& det, size_t j1=0)
  {
    // Same thing, but for Symmetric (not Hermitian) Matrix
    TMVAssert(IsComplex(T()));
    TMVAssert(A.size()>0);
    TMVAssert(!A.isherm());
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(xD.size()+1 == A.size());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);

    const size_t N = A.size();
    static const RealType(T) alpha = (SQRT(RealType(T)(17))+1)/8;

    VectorView<T> D = A.diag();
#ifdef XDEBUG
    Matrix<T> L(N,N,T(0));
    LowerTriMatrixViewOf(L) = A.LowerTri().MakeUnitDiag();
    L.SubMatrix(j1,N,j1,N).SetToIdentity();
    Matrix<T> DD(N,N,T(0));
    DD.diag() = A.diag();
    DD.diag(-1) = xD;
    DD.diag(1) = xD;
    DD.SubMatrix(j1,N,j1,N) = A.SubSymMatrix(j1,N);
    Matrix<T> A0 = L*DD*L.Transpose();
    for(size_t i=j1;i<A.size();i++) P[i]=i;
    A0.ReversePermuteRows(P);
    A0.ReversePermuteCols(P);
#endif

    for (size_t j=j1; j<N;) // ++j or j+=2 done below
    {
      bool seq1 = true;

      if (j==N-1) {
	// No permutation
	P[j] = j;
      } else {
	RealType(T) ajj = abs(D(j));
	size_t p; // p is relative to j index, not absolute.
	RealType(T) apj = A.col(j,j,N).MaxAbsElement(&p);

	if (p == 0 || ajj >= alpha * apj) {
	  // No permutation
	  P[j] = j;
	} else {
	  p+=j;
	  RealType(T) app = abs(D(p));
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
	T dj = D(j);
	if (dj == T(0)) 
	  tmv_error("Zero pivot found in SymLU_Decompose");
	det *= dj;
	A.col(j,j+1,N) /= dj;
	A.SubSymMatrix(j+1,N) -= dj*A.col(j,j+1,N)^A.col(j,j+1,N);
#ifdef XDEBUG
	DD.col(j,j+1,N).Zero();
	DD.row(j,j+1,N).Zero();
	DD(j,j) = D(j);
	L.col(j,j+1,N) = A.col(j,j+1,N);
	DD.SubMatrix(j+1,N,j+1,N) = A.SubSymMatrix(j+1,N);
	Matrix<T> A2 = L*DD*L.Transpose();
	A2.ReversePermuteRows(P);
	A2.ReversePermuteCols(P);
	if (Norm(A0-A2) > 1.e-5*max(RealType(T)(1),Norm(A0))) {
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
	// Invert E:  E^-1 = [ x z ]
	//                   [ z y ]
	T x = D(j);
	T y = D(j+1);
	T z = xD(j) = A(j+1,j);
	A(j+1,j) = T(0);
	T d;
	SymInvert_2x2(x,y,z,&d);
	det *= d;

	if (A.isrm()) {
	  // A(j+2:N,j:j+2) = CE^-1
	  // A(j+2:N,j+2:N) -= CE^-1CT
	  // A(p,q) -= CEinv(p,0)*C(q,0) + CEinv(p,1)*C(q,1)
	  //   p >= q, so loop from bottom up.
	  const int si = A.stepi();
	  T* Cp0 = A.ptr() + (N-1)*si + j;
	  T* Cp1 = Cp0 + 1;
	  for(size_t p = N-1; p>=j+2; --p,Cp0-=si,Cp1-=si) {
	    T Lp0 = *Cp0 * x + *Cp1 * z;
	    T Lp1 = *Cp0 * z + *Cp1 * y;
	    A.row(p,j+2,p+1) -= Lp0 * A.col(j,j+2,p+1);
	    A.row(p,j+2,p+1) -= Lp1 * A.col(j+1,j+2,p+1);
	    *Cp0 = Lp0;
	    *Cp1 = Lp1;
	  }
	} else {
	  // A(j+2:N,j:j+2) = CE^-1
	  // A(j+2:N,j+2:N) -= CE^-1CT 
	  //                -= C (CE^-1)T (since E^-1 is Symmetric)
	  // A(p,q) -= C(p,0)*CEinv(q,0) + C(p,1)*CEinv(q,1)
	  //   p >= q, so loop from top down
	  const int sj = A.stepj();
	  T* Cq0 = A.ptr() + (j+2) + j*sj;
	  T* Cq1 = Cq0 + sj;
	  for(size_t q = j+2; q<N; ++q,++Cq0,++Cq1) {
	    T Lq0 = *Cq0 * x + *Cq1 * z;
	    T Lq1 = *Cq0 * z + *Cq1 * y;
	    A.col(q,q,N) -= A.col(j,q,N) * Lq0;
	    A.col(q,q,N) -= A.col(j+1,q,N) * Lq1;
	    *Cq0 = Lq0;
	    *Cq1 = Lq1;
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
	DD(j,j+1) = xD(j);
	L.col(j,j+2,N) = A.col(j,j+2,N);
	L.col(j+1,j+2,N) = A.col(j+1,j+2,N);
	DD.SubMatrix(j+2,N,j+2,N) = A.SubSymMatrix(j+2,N);
	Matrix<T> A2 = L*DD*L.Transpose();
	A2.ReversePermuteRows(P);
	A2.ReversePermuteCols(P);
	if (Norm(A0-A2) > 1.e-5*max(RealType(T)(1),Norm(A0))) {
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
    DD.diag(1) = xD;
    DD.diag(-1) = xD;
    L = A.LowerTri().MakeUnitDiag();
    Matrix<T> A2 = L*DD*L.Transpose();
    A2.ReversePermuteRows(P);
    A2.ReversePermuteCols(P);

    if (Norm(A2-A0) > 1.e-5*max(RealType(T)(1),SQR(Norm(L))*Norm(DD))) {
      cerr<<"NonBlockSymLUDecompose\n";
      cerr<<"A0 = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"D = "<<D<<endl;
      cerr<<"xD = "<<xD<<endl;
      cerr<<"L = "<<L<<endl;
      cerr<<"DD = "<<DD<<endl;
      cerr<<"L*D*LT = "<<L*DD*L.Transpose()<<endl;
      cerr<<"A2 = "<<A2<<endl;
      cerr<<"A2-A0 = "<<Matrix<T>(A2-A0).Clip(1.e-5)<<endl;
      cerr<<"Norm(A0-A2) = "<<Norm(A0-A2)<<endl;
      cerr<<"nm = "<<Norm(A2-A0)/(SQR(Norm(L))*Norm(DD))<<endl;
      abort();
    }
#endif
  }

  template <class T> void BlockHermLU_Decompose(
      const SymMatrixView<T>& A, const VectorView<T>& xD, 
      size_t* P, RealType(T)& det)
  {
    // Blocked version with level 3 calls

    TMVAssert(A.size()>0);
    TMVAssert(xD.size()+1 == A.size());
    TMVAssert(A.isherm());
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);

#ifdef XDEBUG
    Matrix<T> A0 = A;
#endif

    static const RealType(T) alpha = (SQRT(RealType(T)(17))+1)/8;
    const size_t N = A.size();

    VectorView<RealType(T)> D = A.diag().Real();
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
	  LD.col(jj,j,N) -= A.SubMatrix(j,N,j1,j) * LD.row(j,0,jj).Conjugate();

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
		LD.row(p,0,jj).Conjugate();

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
		// Do important parts of the A.SwapRowsCols(j,p) call,
		// We will overwrite A.col(j), so don't bother swapping into it.
		D(p) = D(j);
		A.row(p,j+1,p) = A.col(j,j+1,p).Conjugate();
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
		  // Do important parts of the A.SwapRowsCols(j+1,p) call,
		  // We will overwrite A.cols(j,j+1), so don't bother 
		  // swapping into them.
		  D(p) = D(j+1);
		  A.row(p,j+2,p) = A.col(j+1,j+2,p).Conjugate();
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
	    RealType(T) dj = D(j);
	    if (dj == T(0))
	      tmv_error("Zero pivot found in HermLU_Decompose");
	    det *= dj;
	    A.col(j,j+1,N) /= dj;
	    ++j;
	  } else {
	    // LD.cols(j,j+1) now hold L.cols(j,j+1) * E 
	    // Invert E:  E^-1 = [ x z* ]
	    //                   [ z y  ]

	    RealType(T) x = D(j) = REAL(LD(j,jj));
	    RealType(T) y = D(j+1) = REAL(LD(j+1,jj+1));
	    T z = xD(j) = LD(j+1,jj);
	    A(j+1,j)=T(0);

	    RealType(T) d;
	    HermInvert_2x2(x,y,z,&d);
	    det *= d;

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
		*Aq1 = *LDq0 * CONJ(z) + *LDq1 * y;
	      }
	    } else {
	      const int si = A.stepi();
	      T* Aq0 = A.ptr() + (j+2)*si + j;
	      T* Aq1 = Aq0 + 1;
	      for(size_t q = j+2; q<N; ++q,Aq0+=si,Aq1+=si,++LDq0,++LDq1) {
		*Aq0 = *LDq0 * x + *LDq1 * z;
		*Aq1 = *LDq0 * CONJ(z) + *LDq1 * y;
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

	// A22 -= L21 * LD21.Adjoint()
	SymMultMM(T(-1),L21,LD21.Adjoint(),1,A22);

      } else NonBlockHermLU_Decompose(A,xD,P,det,j1);

#ifdef XDEBUG
      Matrix<T> L = A.LowerTri();
      L.Cols(j2,N).Zero();
      L.diag().SetAllTo(T(1));

      Matrix<T> DD = DiagMatrixViewOf(A.diag());
      DD.diag(1) = xD.Conjugate();
      DD.diag(-1) = xD;
      DD.SubMatrix(j2,N,j2,N) = A.SubSymMatrix(j2,N);

      Matrix<T> A2 = L * DD * L.Adjoint();
      for(size_t j=j2;j<N;j++) P[j] = j;
      A2.ReversePermuteRows(P);
      A2.ReversePermuteCols(P);

      if (Norm(A2-A0) > 0.000001*max(RealType(T)(1),SQR(Norm(L))*Norm(DD))) {
	cerr<<"BlockHermLUDecompose\n";
	cerr<<"j1 = "<<j1<<", j2 = "<<j2<<", nb = "<<SYM_LU_BLOCKSIZE<<endl;
	HermMatrix<T,Lower,ColMajor> A1(A.size());
	Vector<T> xD1(xD.size(),T(0));
	size_t P1[A.size()];
	RealType(T) d1(1);;
	A1.LowerTri() = LowerTriMatrixViewOf(A0);
	cerr<<"A0 = "<<Type(A)<<A0<<endl;
	cerr<<"A1 = "<<Type(A1)<<A1<<endl;
	NonBlockHermLU_Decompose(A1.View(),xD1.View(),P1,d1);
	cerr<<"L = "<<L<<endl;
	cerr<<"L1 = "<<A1.LowerTri().MakeUnitDiag()<<endl;
	cerr<<"D = "<<A.diag().Real()<<endl;
	cerr<<"xD = "<<xD<<endl;
	cerr<<"D1 = "<<A1.diag().Real()<<endl;
	cerr<<"xD1 = "<<xD1<<endl;
	cerr<<"P = ";
	for(size_t i=0;i<A.size();i++) cerr<<P[i]<<" ";
	cerr<<"\nP1 = ";
	for(size_t i=0;i<A.size();i++) cerr<<P1[i]<<" ";
	cerr<<endl;
	cerr<<"L*D*Lt = "<<L*DD*L.Adjoint()<<endl;
	cerr<<"A2 = "<<A2<<endl;
	cerr<<"A2-A0 = "<<A2-A0<<endl;
	cerr<<"nm = "<<Norm(A2-A0)/(SQR(Norm(L))*Norm(DD))<<endl;

	abort();
      }
#endif
      j1 = j2;
    }
  }

  template <class T> void NonLapHermLU_Decompose(
      const SymMatrixView<T>& A, const VectorView<T>& xD, 
      size_t* P, RealType(T)& det)
  {
    TMVAssert(A.size()>0);
    TMVAssert(xD.size()+1 == A.size());
    TMVAssert(A.isherm());
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);

    if (A.size() > SYM_LU_BLOCKSIZE) 
      BlockHermLU_Decompose(A,xD,P,det);
    else 
      NonBlockHermLU_Decompose(A,xD,P,det);
  }

  template <class T> void BlockSymLU_Decompose(
      const SymMatrixView<T>& A, const VectorView<T>& xD, 
      size_t* P, T& det)
  {
    TMVAssert(IsComplex(T()));
    TMVAssert(A.size()>0);
    TMVAssert(xD.size()+1 == A.size());
    TMVAssert(!A.isherm());
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);

#ifdef XDEBUG
    Matrix<T> A0 = A;
#endif

    static const RealType(T) alpha = (SQRT(RealType(T)(17))+1)/8;
    const size_t N = A.size();

    VectorView<T> D = A.diag();
    Matrix<T,ColMajor> LD(
	N < SYM_LU_BLOCKSIZE ? 0 : N,
	N < SYM_LU_BLOCKSIZE ? 0 : SYM_LU_BLOCKSIZE+1);

    for (size_t j1=0; j1<N; ) {
      size_t j2 = min(j1+SYM_LU_BLOCKSIZE,N);

      if (j2 < N) {
	size_t j;
	for (j=j1; j<j2;) {
	  bool seq1 = true;
	  size_t jj = j-j1;

	  LD.col(jj,j,N) = A.col(j,j,N);
	  LD.col(jj,j,N) -= A.SubMatrix(j,N,j1,j) * LD.row(j,0,jj);

	  if (j == N-1) {
	    P[j] = j;
	  } else {
	    RealType(T) ajj = abs(LD(j,jj));
	    size_t p; 
	    RealType(T) apj = LD.col(jj,j+1,N).MaxAbsElement(&p);

	    if (ajj >= alpha * apj) {
	      P[j] = j;
	    } else {
	      p+=j+1;

	      LD.col(jj+1,j,p) = A.col(p,j,p);
	      LD.col(jj+1,p,N) = A.col(p,p,N);
	      LD.col(jj+1,j,N) -= A.SubMatrix(j,N,j1,j) * LD.row(p,0,jj);

	      RealType(T) app = abs(LD(p,jj+1));
	      RealType(T) apq = LD.col(jj+1,j,p).MaxAbsElement();
	      if (p+1 < N) {
		RealType(T) apq2 = LD.col(jj+1,p+1,N).MaxAbsElement();
		apq = max(apq,apq2);
	      }
	      if (ajj*apq >= alpha * apj * apj) {
		P[j] = j;
	      } else if (app >= alpha * apq) {
		P[j] = p;
		LD.col(jj,j,N) = LD.col(jj+1,j,N);
		D(p) = D(j);
		A.row(p,j+1,p) = A.col(j,j+1,p);
		A.col(p,p+1,N) = A.col(j,p+1,N);
		Swap(A.row(j,0,j),A.row(p,0,j));
		Swap(LD.row(j,0,jj+1),LD.row(p,0,jj+1));
	      } else {
		seq1 = false;
		P[j] = j;
		P[j+1] = p;
		if (p != j+1) {
		  D(p) = D(j+1);
		  A.row(p,j+2,p) = A.col(j+1,j+2,p);
		  A.col(p,p+1,N) = A.col(j+1,p+1,N);
		  Swap(A.row(j+1,0,j),A.row(p,0,j));
		  Swap(LD.row(j+1,0,jj+2),LD.row(p,0,jj+2));
		}
	      }
	    }
	  }

	  if (seq1) {
	    A.col(j,j,N) = LD.col(jj,j,N);
	    T dj = D(j);
	    if (dj == T(0))
	      tmv_error("Zero pivot found in HermLU_Decompose");
	    det *= dj;
	    A.col(j,j+1,N) /= dj;
	    ++j;
	  } else {
	    T x = D(j) = LD(j,jj);
	    T y = D(j+1) = LD(j+1,jj+1);
	    T z = xD(j) = LD(j+1,jj);
	    A(j+1,j)=T(0);

	    T d;
	    SymInvert_2x2(x,y,z,&d);
	    det *= d;

	    const int ldsj = LD.stepj();
	    T* LDq0 = LD.ptr() + (j+2) + jj*ldsj;
	    T* LDq1 = LDq0 + ldsj;
	    if (A.iscm()) {
	      const int sj = A.stepj();
	      T* Aq0 = A.ptr() + (j+2) + j*sj;
	      T* Aq1 = Aq0 + sj;
	      for(size_t q = j+2; q<N; ++q,++Aq0,++Aq1,++LDq0,++LDq1) {
		*Aq0 = *LDq0 * x + *LDq1 * z;
		*Aq1 = *LDq0 * z + *LDq1 * y;
	      }
	    } else {
	      const int si = A.stepi();
	      T* Aq0 = A.ptr() + (j+2)*si + j;
	      T* Aq1 = Aq0 + 1;
	      for(size_t q = j+2; q<N; ++q,Aq0+=si,Aq1+=si,++LDq0,++LDq1) {
		*Aq0 = *LDq0 * x + *LDq1 * z;
		*Aq1 = *LDq0 * z + *LDq1 * y;
	      }
	    }
	    j+=2;
	  }
	}
	j2 = j;

	SymMatrixView<T> A22 = A.SubSymMatrix(j2,N);
	MatrixView<T> L21 = A.SubMatrix(j2,N,j1,j2);
	MatrixView<T> LD21 = LD.SubMatrix(j2,N,0,j2-j1);

	// A22 -= L21 * LD21.Transpose()
	SymMultMM(T(-1),L21,LD21.Transpose(),1,A22);

      } else NonBlockSymLU_Decompose(A,xD,P,det,j1);

#ifdef XDEBUG
      Matrix<T> L = A.LowerTri();
      L.Cols(j2,N).Zero();
      L.diag().SetAllTo(T(1));

      Matrix<T> DD = DiagMatrixViewOf(A.diag());
      DD.diag(1) = DD.diag(-1) = xD;
      DD.SubMatrix(j2,N,j2,N) = A.SubSymMatrix(j2,N);

      Matrix<T> A2 = L * DD * L.Transpose();
      for(size_t j=j2;j<N;j++) P[j] = j;
      A2.ReversePermuteRows(P);
      A2.ReversePermuteCols(P);

      if (Norm(A2-A0) > 0.001*max(RealType(T)(1),SQR(Norm(L))*Norm(DD))) {
	cerr<<"BlockSymLUDecompose\n";
	cerr<<"j1 = "<<j1<<", j2 = "<<j2<<", nb = "<<SYM_LU_BLOCKSIZE<<endl;
	SymMatrix<T,Lower,ColMajor> A1(A.size());
	Vector<T> xD1(xD.size(),T(0));
	size_t P1[A.size()];
	T d1(1);;
	A1.LowerTri() = LowerTriMatrixViewOf(A0);
	cerr<<"A0 = "<<Type(A)<<A0<<endl;
	cerr<<"A1 = "<<Type(A1)<<A1<<endl;
	NonBlockSymLU_Decompose(A1.View(),xD1.View(),P1,d1);
	cerr<<"L = "<<L<<endl;
	cerr<<"L1 = "<<A1.LowerTri().MakeUnitDiag()<<endl;
	cerr<<"D = "<<A.diag().Real()<<endl;
	cerr<<"xD = "<<xD<<endl;
	cerr<<"D1 = "<<A1.diag().Real()<<endl;
	cerr<<"xD1 = "<<xD1<<endl;
	cerr<<"P = ";
	for(size_t i=0;i<A.size();i++) cerr<<P[i]<<" ";
	cerr<<"\nP1 = ";
	for(size_t i=0;i<A.size();i++) cerr<<P1[i]<<" ";
	cerr<<endl;
	cerr<<"L*D*LT = "<<L*DD*L.Transpose()<<endl;
	cerr<<"A2 = "<<A2<<endl;
	cerr<<"nm = "<<Norm(A2-A0)/(SQR(Norm(L))*Norm(DD))<<endl;

	abort();
      }
#endif
      j1 = j2;
    }
  }

  template <class T> void NonLapSymLU_Decompose(
      const SymMatrixView<T>& A, const VectorView<T>& xD, 
      size_t* P, T& det)
  {
    TMVAssert(IsComplex(T()));
    TMVAssert(A.size()>0);
    TMVAssert(xD.size()+1 == A.size());
    TMVAssert(!A.isherm());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);

    if (A.size() > SYM_LU_BLOCKSIZE) 
      BlockSymLU_Decompose(A,xD,P,det);
    else 
      NonBlockSymLU_Decompose(A,xD,P,det);
  }

#ifdef LAP
  template <class T> inline void LapHermLU_Decompose(
      const SymMatrixView<T>& A, const VectorView<T>& xD,
      size_t* P, RealType(T)& det)
  { NonLapHermLU_Decompose(A,xD,P,det); }
  template <class T> inline void LapSymLU_Decompose(
      const SymMatrixView<T>& A, const VectorView<T>& xD,
      size_t* P, T& det)
  { NonLapSymLU_Decompose(A,xD,P,det); }
  template <> inline void LapHermLU_Decompose(
      const SymMatrixView<double>& A, const VectorView<double>& xD,
      size_t* P, double& det)
  {
    TMVAssert(A.size()>0);
    TMVAssert(xD.size()+1 == A.size());
    TMVAssert(A.isherm());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);
    TMVAssert(A.iscm());
    char u = 'L';
    int n = A.size();
    int lda = A.stepj();
    int lap_p[A.size()];
    int lwork = n*LAP_BLOCKSIZE;
    double* work = LAP_DWork(lwork);
    int info;
    dsytrf(&u,&n,A.ptr(),&lda,lap_p,work,&lwork,&info);
    if (info < 0) tmv_error("dsytrf returned info < 0");
    int* pi = lap_p;
    for(size_t i=0;i<A.size();++i,++pi) {
      if (*pi-1 != int(i)) {
	if (*pi > 0) {
	  Swap(A.row(i,0,i),A.row(*pi-1,0,i));
	  P[i] = *pi-1;
	  det *= A(i,i);
	} else {
	  TMVAssert(*pi < 0);
	  TMVAssert(i+1 < A.size());
	  TMVAssert(*(pi+1) == *pi);
	  P[i] = i;
	  P[i+1] = -*pi-1;
	  if (int(i)+1 != -*pi-1) {
	    Swap(A.row(i+1,0,i),A.row(-*pi-1,0,i));
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
  template <> inline void LapHermLU_Decompose(
      const SymMatrixView<complex<double> >& A,
      const VectorView<complex<double> >& xD,
      size_t* P, double& det)
  {
    TMVAssert(A.size()>0);
    TMVAssert(xD.size()+1 == A.size());
    TMVAssert(A.isherm());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);
    TMVAssert(A.iscm());
    TMVAssert(A.isherm());
    char u = 'L';
    int n = A.size();
    int lda = A.stepj();
    int lap_p[A.size()];
    int lwork = n*LAP_BLOCKSIZE;
    complex<double>* work = LAP_ZWork(lwork);
    int info;
    zhetrf(&u,&n,LAP_Complex(A.ptr()),&lda,lap_p,LAP_Complex(work),&lwork,
	&info);
    if (info < 0) tmv_error("zhetrf returned info < 0");
    int* pi = lap_p;
    for(size_t i=0;i<A.size();++i,++pi) {
      if (*pi-1 != int(i)) {
	if (*pi > 0) {
	  Swap(A.row(i,0,i),A.row(*pi-1,0,i));
	  P[i] = *pi-1;
	  det *= REAL(A(i,i));
	} else {
	  TMVAssert(*pi < 0);
	  TMVAssert(i+1 < A.size());
	  TMVAssert(*(pi+1) == *pi);
	  P[i] = i;
	  P[i+1] = -*pi-1;
	  if (int(i)+1 != -*pi-1) {
	    Swap(A.row(i+1,0,i),A.row(-*pi-1,0,i));
	  }
	  complex<double>& xDi = A(i+1,i).GetRef();
	  xD(i) = xDi;
	  det *= REAL(A(i,i))*REAL(A(i+1,i+1))-NORM(xDi);
	  xDi = double(0);
	  ++i; ++pi; // extra ++ for 2x2 case
	}
      } else {
	det *= REAL(A(i,i));
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
    TMVAssert(!A.isherm());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);
    TMVAssert(A.iscm());
    TMVAssert(!A.isherm());
    char u = 'L';
    int n = A.size();
    int lda = A.stepj();
    int lap_p[A.size()];
    int lwork = n*LAP_BLOCKSIZE;
    complex<double>* work = LAP_ZWork(lwork);
    int info;
    zsytrf(&u,&n,LAP_Complex(A.ptr()),&lda,lap_p,LAP_Complex(work),&lwork,
	&info);
    if (info < 0) tmv_error("zsytrf returned info < 0");
    int* pi = lap_p;
    for(size_t i=0;i<A.size();++i,++pi) {
      if (*pi-1 != int(i)) {
	if (*pi > 0) {
	  Swap(A.row(i,0,i),A.row(*pi-1,0,i));
	  P[i] = *pi-1;
	  det *= A(i,i);
	} else {
	  TMVAssert(*pi < 0);
	  TMVAssert(i+1 < A.size());
	  TMVAssert(*(pi+1) == *pi);
	  P[i] = i;
	  P[i+1] = -*pi-1;
	  if (int(i)+1 != -*pi-1) {
	    Swap(A.row(i+1,0,i),A.row(-*pi-1,0,i));
	  }
	  complex<double>& xDi = A(i+1,i).GetRef();
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
#ifndef NOFLOAT
  template <> inline void LapHermLU_Decompose(
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
    char u = 'L';
    int n = A.size();
    int lda = A.stepj();
    int lap_p[A.size()];
    int lwork = n*LAP_BLOCKSIZE;
    float* work = LAP_SWork(lwork);
    int info;
    ssytrf(&u,&n,A.ptr(),&lda,lap_p,work,&lwork,&info);
    if (info < 0) tmv_error("ssytrf returned info < 0");
    int* pi = lap_p;
    for(size_t i=0;i<A.size();++i,++pi) {
      if (*pi-1 != int(i)) {
	if (*pi > 0) {
	  Swap(A.row(i,0,i),A.row(*pi-1,0,i));
	  P[i] = *pi-1;
	  det *= A(i,i);
	} else {
	  TMVAssert(*pi < 0);
	  TMVAssert(i+1 < A.size());
	  TMVAssert(*(pi+1) == *pi);
	  P[i] = i;
	  P[i+1] = -*pi-1;
	  if (int(i)+1 != -*pi-1) {
	    Swap(A.row(i+1,0,i),A.row(-*pi-1,0,i));
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
  template <> inline void LapHermLU_Decompose(
      const SymMatrixView<complex<float> >& A,
      const VectorView<complex<float> >& xD,
      size_t* P, float& det)
  {
    TMVAssert(A.size()>0);
    TMVAssert(xD.size()+1 == A.size());
    TMVAssert(A.isherm());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);
    TMVAssert(A.iscm());
    TMVAssert(A.isherm());
    char u = 'L';
    int n = A.size();
    int lda = A.stepj();
    int lap_p[A.size()];
    int lwork = n*LAP_BLOCKSIZE;
    complex<float>* work = LAP_CWork(lwork);
    int info;
    chetrf(&u,&n,LAP_Complex(A.ptr()),&lda,lap_p,LAP_Complex(work),&lwork,
	&info);
    if (info < 0) tmv_error("chetrf returned info < 0");
    int* pi = lap_p;
    for(size_t i=0;i<A.size();++i,++pi) {
      if (*pi-1 != int(i)) {
	if (*pi > 0) {
	  Swap(A.row(i,0,i),A.row(*pi-1,0,i));
	  P[i] = *pi-1;
	  det *= REAL(A(i,i));
	} else {
	  TMVAssert(*pi < 0);
	  TMVAssert(i+1 < A.size());
	  TMVAssert(*(pi+1) == *pi);
	  P[i] = i;
	  P[i+1] = -*pi-1;
	  if (int(i)+1 != -*pi-1) {
	    Swap(A.row(i+1,0,i),A.row(-*pi-1,0,i));
	  }
	  complex<float>& xDi = A(i+1,i).GetRef();
	  xD(i) = xDi;
	  det *= REAL(A(i,i))*REAL(A(i+1,i+1))-NORM(xDi);
	  xDi = float(0);
	  ++i; ++pi; // extra ++ for 2x2 case
	}
      } else {
	det *= REAL(A(i,i));
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
    TMVAssert(!A.isherm());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);
    TMVAssert(A.iscm());
    TMVAssert(!A.isherm());
    char u = 'L';
    int n = A.size();
    int lda = A.stepj();
    int lap_p[A.size()];
    int lwork = n*LAP_BLOCKSIZE;
    complex<float>* work = LAP_CWork(lwork);
    int info;
    csytrf(&u,&n,LAP_Complex(A.ptr()),&lda,lap_p,LAP_Complex(work),&lwork,
	&info);
    if (info < 0) tmv_error("csytrf returned info < 0");
    int* pi = lap_p;
    for(size_t i=0;i<A.size();++i,++pi) {
      if (*pi-1 != int(i)) {
	if (*pi > 0) {
	  Swap(A.row(i,0,i),A.row(*pi-1,0,i));
	  P[i] = *pi-1;
	  det *= A(i,i);
	} else {
	  TMVAssert(*pi < 0);
	  TMVAssert(i+1 < A.size());
	  TMVAssert(*(pi+1) == *pi);
	  P[i] = i;
	  P[i+1] = -*pi-1;
	  if (int(i)+1 != -*pi-1) {
	    Swap(A.row(i+1,0,i),A.row(-*pi-1,0,i));
	  }
	  complex<float>& xDi = A(i+1,i).GetRef();
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
#endif // NOFLOAT
#endif // LAP
  template <class T> inline void SymLU_Decompose(
      const SymMatrixView<T>& A, const VectorView<T>& xD, 
      size_t* P, T& det)
  {
    TMVAssert(A.size() > 0);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(xD.size()+1 == A.size());
    if (A.isconj()) SymLU_Decompose(A.Conjugate(),xD,P,det);
    else if (A.uplo() == Upper) {
      if (A.isherm()) SymLU_Decompose(A.Adjoint(),xD,P,det);
      else SymLU_Decompose(A.Transpose(),xD,P,det);
    }
    else {
      xD.Zero();
#ifdef LAP
      if (A.iscm())
	if (A.isherm()) {
	  RealType(T) det2(1);
	  LapHermLU_Decompose(A,xD,P,det2);
	  det *= det2;
	}
	else LapSymLU_Decompose(A,xD,P,det);
      else
#endif
	if (A.isherm()) {
	  RealType(T) det2(1);
	  NonLapHermLU_Decompose(A,xD,P,det2);
	  det *= det2;
	}
	else NonLapSymLU_Decompose(A,xD,P,det);
    }
  }

#define InstFile "TMV_SymLUDiv_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


