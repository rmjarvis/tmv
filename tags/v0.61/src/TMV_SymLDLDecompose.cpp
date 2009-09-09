///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2008                                                        //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file LICENSE.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



#include "TMV_Blas.h"
#include "TMV_SymLDLDiv.h"
#include "TMV_SymLDLD.h"
#include "TMV_SymMatrix.h"
#include "TMV_SymBandMatrix.h"
#include "TMV_Vector.h"
#include "TMV_VectorArith.h"
#include "TMV_Matrix.h"
#include "TMV_SymMatrixArith.h"
#include "TMV_MatrixArith.h"
#include "TMV_VIt.h"

//#define XDEBUG

#ifdef XDEBUG
#include "TMV_TriMatrixArith.h"
#include "TMV_SymBandMatrixArith.h"
#include "TMV_VIt.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

#define RT RealType(T)

#ifdef TMV_BLOCKSIZE
#define SYM_LDL_BLOCKSIZE TMV_BLOCKSIZE
#else
#define SYM_LDL_BLOCKSIZE 48
#endif

  //
  // Decompose
  //

  // MJ: This is a bit slower for row major, even with blocking, so write the
  // row major version where A is decomposed into Lt D L rather than L D Lt.
  template <bool herm, class T> static void NonBlockLDL_Decompose(
      const SymMatrixView<T>& A, const VectorView<T>& xD, 
      int* P, RT& logdet, T& signdet, int j1=0)
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
    const int N = A.size();
    static const RT alpha = (SQRT(RT(17))+1)/8;

#ifdef XTEST 
    TMVAssert(A.HermOK());
#endif

    VectorView<T> D = A.diag();
    T* Dj = D.ptr()+j1*D.step();
    for (int j=j1; j<N;) { // ++j or j+=2 done below
      //cout<<"j = "<<j<<endl;
      //cout<<"A = "<<A<<endl;
      bool seq1 = true;
      if (j == N-1) {
	// No permutation
	P[j] = j;
      } else {
	RT ajj = herm ? ABS(REAL(*Dj)) : ABS(*Dj);
	int p; // p is relative to j index, not absolute.
	RT apj = A.col(j,j,N).MaxAbsElement(&p);

	if (p == 0 || ajj >= alpha * apj) {
	  // No permutation
	  P[j] = j;
	} else {
	  p+=j;
	  T* Dp = D.ptr()+p*D.step();
	  RT app = herm ? ABS(REAL(*Dp)) : ABS(*Dp);
	  RT apq = A.row(p,j,p).MaxAbsElement();
	  if (p+1 < N) {
	    RT apq2 = A.col(p,p+1,N).MaxAbsElement();
	    apq = MAX(apq,apq2);
	  }
	  if (ajj*apq >= alpha * apj * apj) {
	    // No permutation
	    P[j] = j;
	  } else if (app >= alpha * apq) {
	    // Permute p diagonal into j spot
	    TMVAssert(p<int(A.size()));
	    A.SwapRowsCols(j,p);
	    P[j] = p;
	  } else {
	    // Permute pj element into j+1,j spot
	    // This also permutes pp element into j+1,j+1
	    seq1 = false;
	    P[j] = j;
	    P[j+1] = p;
	    if (p != j+1) {
	      TMVAssert(p<int(A.size()));
	      A.SwapRowsCols(j+1,p);
	    }
	  }
	}
      }

      //cout<<"Before LDL solving: A = "<<A<<endl;
      // Now the LDL solving:
      if (seq1) {
	if (herm) {
	  RT dj = REAL(*Dj);
	  if (signdet != T(0)) {
	    if (dj == RT(0)) signdet = T(0);
	    else {
	      if (dj < RT(0)) signdet = -signdet;
	      logdet += std::log(ABS(dj));
	    }
	  }
	  if (dj != RT(0)) {
	    //cout<<"0\n";
	    //cout<<"A = "<<A<<endl;
	    A.col(j,j+1,N) /= dj;
	    //cout<<"1\n";
	    //cout<<"A = "<<A<<endl;
	    A.SubSymMatrix(j+1,N) -= 
	      dj*A.col(j,j+1,N)^A.col(j,j+1,N).Conjugate();
	    //cout<<"2\n";
	    //cout<<"A = "<<A<<endl;
	  }
	} else {
	  T dj = *Dj;
	  if (signdet != T(0)) {
	    if (dj == T(0)) signdet = T(0);
	    else {
	      RT absd = ABS(dj);
	      signdet *= dj/absd;
	      logdet += std::log(absd);
	    }
	  }
	  if (dj != T(0)) {
	    //cout<<"0b\n";
	    //cout<<"A = "<<A<<endl;
	    A.col(j,j+1,N) /= dj;
	    //cout<<"3\n";
	    //cout<<"A = "<<A<<endl;
	    A.SubSymMatrix(j+1,N) -= dj*A.col(j,j+1,N)^A.col(j,j+1,N);
	    //cout<<"4\n";
	    //cout<<"A = "<<A<<endl;
	  }
	}
	++j; Dj += D.step();
      } else {
	// Invert E:  E^-1 = [ x z* ]
	//                   [ z y  ]

	T x = *Dj;
	// Dj is at A(j,j), so
	// A(j+1,j) = *(Dj+A.stepi())
	T* Ax = Dj + A.stepi();
	T z = xD[j] = *Ax;
	T y = *(Dj+=D.step());
#ifdef TMVFLDEBUG
	TMVAssert(Ax >= A.first);
	TMVAssert(Ax< A.last);
#endif
	*Ax=T(0);
	T d;
	SymInvert_2x2<herm>(x,y,z,&d);
	if (signdet != T(0)) {
	  if (d == T(0)) signdet = T(0);
	  else if (herm) {
	    if (REAL(d) < RT(0)) signdet = -signdet; 
	    logdet += std::log(ABS(REAL(d)));
	  } else {
	    RT absd = ABS(d);
	    signdet *= d/absd;
	    logdet += std::log(absd);
	  }
	}
#ifdef XTEST
	RT NormEinv = NORM(x) + NORM(y) + 2*NORM(z);
	if (NormEinv < 1.) NormEinv = 1.;
	RT NormC = Norm(A.SubMatrix(j+2,N,j,j+2));
	if (NormC < 1.) NormC = 1.;
	RT eps = NormC*NormC*NormEinv;
	eps *= N*Epsilon<T>();
#endif

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
	  for(int q = j+2; q<N; ++q,++Cq0,++Cq1) {
	    T Lq0 = *Cq0 * x + *Cq1 * z;
	    T Lq1 = *Cq0 * (herm?CONJ(z):z) + *Cq1 * y;
	    A.col(q,q,N) -= A.col(j,q,N) * (herm?CONJ(Lq0):Lq0);
	    A.col(q,q,N) -= A.col(j+1,q,N) * (herm?CONJ(Lq1):Lq1);
#ifdef TMVFLDEBUG
	    TMVAssert(Cq0 >= A.first);
	    TMVAssert(Cq0< A.last);
	    TMVAssert(Cq1 >= A.first);
	    TMVAssert(Cq1< A.last);
#endif
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
	  for(int p = N-1; p>=j+2; --p,Cp0-=si,Cp1-=si) {
	    T Lp0 = *Cp0 * x + *Cp1 * z;
	    T Lp1 = *Cp0 * (herm?CONJ(z):z) + *Cp1 * y;
	    A.row(p,j+2,p+1) -= Lp0 * 
	      (herm ? A.col(j,j+2,p+1).Conjugate() : A.col(j,j+2,p+1));
	    A.row(p,j+2,p+1) -= Lp1 *
	      (herm ? A.col(j+1,j+2,p+1).Conjugate() : A.col(j+1,j+2,p+1));
#ifdef TMVFLDEBUG
	    TMVAssert(Cp0 >= A.first);
	    TMVAssert(Cp0< A.last);
	    TMVAssert(Cp1 >= A.first);
	    TMVAssert(Cp1< A.last);
#endif
	    *Cp0 = Lp0;
	    *Cp1 = Lp1;
	  }
	}
	if (herm && IsComplex(T())) {
#ifdef XTEST
	  TMVAssert(NormInf(A.diag().Imag()) <= eps);
#endif
	  A.diag().SubVector(j+2,N).Imag().Zero();
	}

	j+=2; Dj += D.step(); // Already did one += D.step() above
      }
    }
#ifdef XTEST 
    TMVAssert(A.HermOK());
#endif
  }

  template <bool herm, class T> static void BlockLDL_Decompose(
      const SymMatrixView<T>& A, const VectorView<T>& xD, 
      int* P, RT& logdet, T& signdet)
  {
    // Blocked version with level 3 calls

    TMVAssert(A.size()>0);
    TMVAssert(xD.size()+1 == A.size());
    TMVAssert(IsReal(T()) || herm == A.isherm());
    TMVAssert(herm || IsComplex(T()));
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);

#ifdef XTEST
    TMVAssert(A.HermOK());
#endif

    static const RT alpha = (SQRT(RT(17))+1)/8;
    const int N = A.size();

    VectorView<T> D = A.diag();
    // Only make LD if we're going to use it:
    Matrix<T,ColMajor> LD(
	N < SYM_LDL_BLOCKSIZE ? 0 : N,
	N < SYM_LDL_BLOCKSIZE ? 0 : SYM_LDL_BLOCKSIZE+1);
    T* Dj = D.ptr();
    for (int j1=0; j1<N; ) {
      int j2 = MIN(j1+SYM_LDL_BLOCKSIZE,N);

      if (j2 < N) { // On last loop we skip some steps.  See below.
	int j;
	T* LDjj = LD.ptr()+j1;  // = LD(j,jj)
	for (j=j1; j<j2;) { // ++j or j+=2 done below
	  bool seq1 = true;
	  int jj = j-j1; // Col index in LD

	  LD.col(jj,j,N) = A.col(j,j,N);
	  LD.col(jj,j,N) -= A.SubMatrix(j,N,j1,j) *
	    (herm ? LD.row(j,0,jj).Conjugate() : LD.row(j,0,jj));

	  if (j == N-1) {
	    // No permutation
	    P[j] = j;
	  } else {
	    RT ajj = ABS(*LDjj);
	    int p; // p is relative to j+1 index, not absolute.
	    RT apj = LD.col(jj,j+1,N).MaxAbsElement(&p);

	    if (ajj >= alpha * apj) {
	      // No permutation
	      P[j] = j;
	    } else {
	      p+=j+1;

	      LD.col(jj+1,j,p) = A.col(p,j,p);
	      LD.col(jj+1,p,N) = A.col(p,p,N);
	      LD.col(jj+1,j,N) -= A.SubMatrix(j,N,j1,j) * 
		(herm ? LD.row(p,0,jj).Conjugate() : LD.row(p,0,jj));

	      //RT app = ABS(LD(p,jj+1));
	      const T* LDpj = LD.cptr() + p + (j1+1)*LD.stepj();
	      RT app = ABS(*LDpj);
	      RT apq = LD.col(jj+1,j,p).MaxAbsElement();
	      if (p+1 < N) {
		RT apq2 = LD.col(jj+1,p+1,N).MaxAbsElement();
		apq = MAX(apq,apq2);
	      }
	      if (ajj*apq >= alpha * apj * apj) {
		// No permutation
		P[j] = j;
	      } else if (app >= alpha * apq) {
		// Permute p diagonal into j spot
		P[j] = p;
		// Move updated pivot column (in LD) to j column
		LD.col(jj,j,N) = LD.col(jj+1,j,N);
		TMVAssert(p<int(A.size()));
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
		  TMVAssert(p<int(A.size()));
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

	  // Now the LDL solving:
	  if (seq1) {
	    // LD.col(j) now holds L.col(j)*D(j)
	    A.col(j,j,N) = LD.col(jj,j,N);
	    if (herm) {
	      if (IsComplex(T())) {
#ifdef TMVFLDEBUG
		TMVAssert(D.ptr() >= A.first);
		TMVAssert(D.ptr() < A.last);
#endif
		*Dj = REAL(*Dj);
	      }

	      RT dj = REAL(*Dj);
	      if (signdet != T(0)) {
		if (dj == RT(0)) signdet = T(0);
		else {
		  if (dj < RT(0)) signdet = -signdet;
		  logdet += std::log(ABS(dj));
		}
	      }
	      if (dj != RT(0)) A.col(j,j+1,N) /= dj;
	    } else {
	      T dj = *Dj;
	      if (signdet != T(0)) {
		if (dj == RT(0)) signdet = T(0);
		else {
		  RT ad = ABS(dj);
		  signdet *= dj/ad;
		  logdet += std::log(ad);
		}
	      }
	      if (dj != T(0)) A.col(j,j+1,N) /= dj;
	    }
	    ++j; 
	    Dj+=D.step();
	    LDjj += LD.stepj()+1;
	  } else {
	    // LD.cols(j,j+1) now hold L.cols(j,j+1) * E 
	    // Invert E:  E^-1 = [ x z* ]
	    //                   [ z y  ]

#ifdef TMVFLDEBUG
	    TMVAssert(Dj >= A.first);
	    TMVAssert(Dj < A.last);
#endif
	    if (herm) *Dj = REAL(*LDjj);
	    else *Dj = *LDjj;
	    T x = *Dj;

	    ++LDjj; // Now LD(j+1,jj+1)
	    T z = xD[j] = *LDjj;

	    Dj += D.step(); // Now D(j+1)
	    LDjj += LD.stepj(); // Now LD(j+1,jj+1)
#ifdef TMVFLDEBUG
	    TMVAssert(Dj >= A.first);
	    TMVAssert(Dj < A.last);
#endif
	    if (herm) *Dj = REAL(*LDjj);
	    else *Dj = *LDjj;
	    T y = *Dj; // = D(j+1)

	    A(j+1,j)=T(0);

	    T d;
	    SymInvert_2x2<herm>(x,y,z,&d);
	    if (signdet != T(0)) {
	      if (d == RT(0)) signdet = T(0);
	      else if (herm) {
		RT ad = ABS(REAL(d));
		if (REAL(d) < 0) signdet = -signdet;
		logdet += std::log(ad);
	      } else {
		RT ad = ABS(d);
		signdet *= d/ad;
		logdet += std::log(ad);
	      }
	    }

	    // A(j+2:N,j:j+2) = LD E^-1
	    const int ldsj = LD.stepj();
	    T* LDq0 = LD.ptr() + (j+2) + jj*ldsj;
	    T* LDq1 = LDq0 + ldsj;
	    if (A.iscm()) {
	      const int sj = A.stepj();
	      T* Aq0 = A.ptr() + (j+2) + j*sj;
	      T* Aq1 = Aq0 + sj;
	      for(int q = j+2; q<N; ++q,++Aq0,++Aq1,++LDq0,++LDq1) {
#ifdef TMVFLDEBUG
		TMVAssert(Aq0 >= A.first);
		TMVAssert(Aq0 < A.last);
		TMVAssert(Aq1 >= A.first);
		TMVAssert(Aq1 < A.last);
#endif
		*Aq0 = *LDq0 * x + *LDq1 * z;
		*Aq1 = *LDq0 * (herm?CONJ(z):z) + *LDq1 * y;
	      }
	    } else {
	      const int si = A.stepi();
	      T* Aq0 = A.ptr() + (j+2)*si + j;
	      T* Aq1 = Aq0 + 1;
	      for(int q = j+2; q<N; ++q,Aq0+=si,Aq1+=si,++LDq0,++LDq1) {
#ifdef TMVFLDEBUG
		TMVAssert(Aq0 >= A.first);
		TMVAssert(Aq0 < A.last);
		TMVAssert(Aq1 >= A.first);
		TMVAssert(Aq1 < A.last);
#endif
		*Aq0 = *LDq0 * x + *LDq1 * z;
		*Aq1 = *LDq0 * (herm?CONJ(z):z) + *LDq1 * y;
	      }
	    }
	    j+=2; 
	    Dj+=D.step(); // One of these steps already done above
	    LDjj += LD.stepj()+1;
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
	SymMultMM<true>(T(-1),L21,herm?LD21.Adjoint():LD21.Transpose(),A22);

      } else NonBlockLDL_Decompose<herm>(A,xD,P,logdet,signdet,j1);

      j1 = j2;
    }
#ifdef XTEST
    TMVAssert(A.HermOK());
#endif
  }

  template <bool herm, class T> static void NonLapLDL_Decompose(
      const SymMatrixView<T>& A, const VectorView<T>& xD, 
      int* P, RT& logdet, T& signdet)
  {
    TMVAssert(A.size()>0);
    TMVAssert(xD.size()+1 == A.size());
    TMVAssert(herm == A.isherm());
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);

    if (A.size() > SYM_LDL_BLOCKSIZE) 
      BlockLDL_Decompose<herm>(A,xD,P,logdet,signdet);
    else 
      NonBlockLDL_Decompose<herm>(A,xD,P,logdet,signdet);
  }

#ifdef LAP
  template <class T> static inline void LapLDL_Decompose(
      const SymMatrixView<T>& A, const VectorView<T>& xD,
      int* P, RT& logdet, T& signdet)
  { 
    if (A.isherm())
      NonLapLDL_Decompose<true>(A,xD,P,logdet,signdet); 
    else
      NonLapLDL_Decompose<false>(A,xD,P,logdet,signdet); 
  }
#ifdef INST_DOUBLE
  template <> void LapLDL_Decompose(
      const SymMatrixView<double>& A, const VectorView<double>& xD,
      int* P, double& logdet, double& signdet)
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
    double* Aii = A.ptr();
    const int Astep = A.stepj()+1;
    for(int i=0;i<n;++i,++pi,Aii+=Astep) {
      int ip = *pi LAPMINUS1;
      if (ip != int(i)) {
	if (ip >= 0) {
	  Swap(A.row(i,0,i),A.row(ip,0,i));
	  P[i] = ip;
	  if (signdet != double(0)) {
	    if (*Aii == double(0)) signdet = double(0);
	    else {
	      if (*Aii < double(0)) signdet = -signdet;
	      logdet += std::log(ABS(*Aii));
	    }
	  }
	} else {
	  TMVAssert(i+1 < A.size());
	  TMVAssert(*(pi+1) == *pi);
	  ip = -*pi LAPMINUS1;
	  P[i] = i;
	  P[i+1] = ip;
	  if (int(i+1) != ip) {
	    Swap(A.row(i+1,0,i),A.row(ip,0,i));
	  }
	  double a = *Aii;
	  double* Aoff = Aii+1;
	  double b = *Aoff;
	  *Aoff = double(0);
	  Aii += Astep;
	  double c = *Aii;

	  xD[i] = b;
	  double d = a*c-b*b;
	  if (signdet != double(0)) {
	    if (d == double(0)) signdet = double(0);
	    else {
	      if (d < double(0)) signdet = -signdet;
	      logdet += std::log(ABS(d));
	    }
	  }
	  ++i; ++pi; // extra ++ for 2x2 case
	}
      } else {
	if (signdet != double(0)) {
	  if (*Aii == double(0)) signdet = double(0);
	  else {
	    if (*Aii < double(0)) signdet = -signdet;
	    logdet += std::log(ABS(*Aii));
	  }
	}
	P[i] = i;
      }
    }
  }
  template <> void LapLDL_Decompose(
      const SymMatrixView<std::complex<double> >& A,
      const VectorView<std::complex<double> >& xD,
      int* P, double& logdet, std::complex<double>& signdet)
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
    std::complex<double>* work = LAP_ZWork(lwork);
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
    std::complex<double>* Aii = A.ptr();
    const int Astep = A.stepj()+1;
    for(int i=0;i<n;++i,++pi,Aii+=Astep) {
      int ip = *pi LAPMINUS1;
      if (ip != int(i)) {
	if (ip >= 0) {
	  Swap(A.row(i,0,i),A.row(ip,0,i));
	  P[i] = ip;
	  if (signdet != double(0)) {
	    if (*Aii == double(0)) signdet = double(0);
	    else if (A.isherm()) {
	      if (REAL(*Aii) < double(0)) signdet = -signdet;
	      logdet += std::log(ABS(REAL(*Aii)));
	    } else {
	      double ad = ABS(*Aii);
	      signdet *= *Aii/ad;
	      logdet += std::log(ad);
	    }
	  }
	} else {
	  TMVAssert(i+1 < A.size());
	  TMVAssert(*(pi+1) == *pi);
	  ip = -*pi LAPMINUS1;
	  P[i] = i;
	  P[i+1] = ip;
	  if (int(i+1) != ip) {
	    Swap(A.row(i+1,0,i),A.row(ip,0,i));
	  }
	  std::complex<double>* Aoff = Aii+1;
	  std::complex<double> b = *Aoff;
	  *Aoff = double(0);
	  xD[i] = b;
	  if (signdet != double(0)) {
	    if (A.isherm()) {
	      double a = REAL(*Aii);
	      Aii += Astep;
	      double c = REAL(*Aii);
	      double d = a*c-NORM(b);
	      if (d == double(0)) signdet = double(0);
	      else {
		if (d < double(0)) signdet = -signdet;
		logdet += std::log(ABS(d));
	      }
	    } else {
	      std::complex<double> a = *Aii;
	      Aii += Astep;
	      std::complex<double> c = *Aii;
	      std::complex<double> d = a*c-b*b;
	      if (d == double(0)) signdet = double(0);
	      else {
		double ad = ABS(d);
		signdet *= d/ad;
		logdet += std::log(ad);
	      }
	    }
	  } else { Aii += Astep; }
	  ++i; ++pi; // extra ++ for 2x2 case
	}
      } else {
	if (signdet != double(0)) {
	  if (*Aii == double(0)) signdet = double(0);
	  else if (A.isherm()) {
	    if (REAL(*Aii) < double(0)) signdet = -signdet;
	    logdet += std::log(ABS(REAL(*Aii)));
	  } else {
	    double ad = ABS(*Aii);
	    signdet *= *Aii/ad;
	    logdet += std::log(ad);
	  }
	}
	P[i] = i;
      }
    }
  }
#endif
#ifdef INST_FLOAT
  template <> void LapLDL_Decompose(
      const SymMatrixView<float>& A, const VectorView<float>& xD,
      int* P, float& logdet, float& signdet)
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
    float* work = LAP_SWork(lwork);
#endif
    LAPNAME(ssytrf) (LAPCM LAPCH_LO,LAPV(n),LAPP(A.ptr()),LAPV(lda),
	LAPP(lap_p.get()) LAPWK(work) LAPVWK(lwork) LAPINFO LAP1);
    LAP_Results("ssytrf");
    int* pi = lap_p.get();
    float* Aii = A.ptr();
    const int Astep = A.stepj()+1;
    for(int i=0;i<n;++i,++pi,Aii+=Astep) {
      int ip = *pi LAPMINUS1;
      if (ip != int(i)) {
	if (ip >= 0) {
	  Swap(A.row(i,0,i),A.row(ip,0,i));
	  P[i] = ip;
	  if (signdet != float(0)) {
	    if (*Aii == float(0)) signdet = float(0);
	    else {
	      if (*Aii < float(0)) signdet = -signdet;
	      logdet += std::log(ABS(*Aii));
	    }
	  }
	} else {
	  TMVAssert(i+1 < A.size());
	  TMVAssert(*(pi+1) == *pi);
	  ip = -*pi LAPMINUS1;
	  P[i] = i;
	  P[i+1] = ip;
	  if (int(i+1) != ip) {
	    Swap(A.row(i+1,0,i),A.row(ip,0,i));
	  }
	  float a = *Aii;
	  float* Aoff = Aii+1;
	  float b = *Aoff;
	  *Aoff = float(0);
	  Aii += Astep;
	  float c = *Aii;

	  xD[i] = b;
	  float d = a*c-b*b;
	  if (signdet != float(0)) {
	    if (d == float(0)) signdet = float(0);
	    else {
	      if (d < float(0)) signdet = -signdet;
	      logdet += std::log(ABS(d));
	    }
	  }
	  ++i; ++pi; // extra ++ for 2x2 case
	}
      } else {
	if (signdet != float(0)) {
	  if (*Aii == float(0)) signdet = float(0);
	  else {
	    if (*Aii < float(0)) signdet = -signdet;
	    logdet += std::log(ABS(*Aii));
	  }
	}
	P[i] = i;
      }
    }
  }
  template <> void LapLDL_Decompose(
      const SymMatrixView<std::complex<float> >& A,
      const VectorView<std::complex<float> >& xD,
      int* P, float& logdet, std::complex<float>& signdet)
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
    std::complex<float>* work = LAP_CWork(lwork);
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
    std::complex<float>* Aii = A.ptr();
    const int Astep = A.stepj()+1;
    for(int i=0;i<n;++i,++pi,Aii+=Astep) {
      int ip = *pi LAPMINUS1;
      if (ip != int(i)) {
	if (ip >= 0) {
	  Swap(A.row(i,0,i),A.row(ip,0,i));
	  P[i] = ip;
	  if (signdet != float(0)) {
	    if (*Aii == float(0)) signdet = float(0);
	    else if (A.isherm()) {
	      if (REAL(*Aii) < float(0)) signdet = -signdet;
	      logdet += std::log(ABS(REAL(*Aii)));
	    } else {
	      float ad = ABS(*Aii);
	      signdet *= *Aii/ad;
	      logdet += std::log(ad);
	    }
	  }
	} else {
	  TMVAssert(i+1 < A.size());
	  TMVAssert(*(pi+1) == *pi);
	  ip = -*pi LAPMINUS1;
	  P[i] = i;
	  P[i+1] = ip;
	  if (int(i+1) != ip) {
	    Swap(A.row(i+1,0,i),A.row(ip,0,i));
	  }
	  std::complex<float>* Aoff = Aii+1;
	  std::complex<float> b = *Aoff;
	  *Aoff = float(0);
	  xD[i] = b;
	  if (signdet != float(0)) {
	    if (A.isherm()) {
	      float a = REAL(*Aii);
	      Aii += Astep;
	      float c = REAL(*Aii);
	      float d = a*c-NORM(b);
	      if (d == float(0)) signdet = float(0);
	      else {
		if (d < float(0)) signdet = -signdet;
		logdet += std::log(ABS(d));
	      }
	    } else {
	      std::complex<float> a = *Aii;
	      Aii += Astep;
	      std::complex<float> c = *Aii;
	      std::complex<float> d = a*c-b*b;
	      if (d == float(0)) signdet = float(0);
	      else {
		float ad = ABS(d);
		signdet *= d/ad;
		logdet += std::log(ad);
	      }
	    }
	  } else { Aii += Astep; }
	  ++i; ++pi; // extra ++ for 2x2 case
	}
      } else {
	if (signdet != float(0)) {
	  if (*Aii == float(0)) signdet = float(0);
	  else if (A.isherm()) {
	    if (REAL(*Aii) < float(0)) signdet = -signdet;
	    logdet += std::log(ABS(REAL(*Aii)));
	  } else {
	    float ad = ABS(*Aii);
	    signdet *= *Aii/ad;
	    logdet += std::log(ad);
	  }
	}
	P[i] = i;
      }
    }
  }
#endif
#endif // LAP

  template <class T> void LDL_Decompose(
      const SymMatrixView<T>& A, const VectorView<T>& xD, 
      int* P, RT& logdet, T& signdet)
  {
    TMVAssert(A.size() > 0);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(xD.size()+1 == A.size());
#ifdef XDEBUG
    Matrix<T> A0(A);
    //cout<<"Start LDL_Decomp\n";
#endif
#ifdef XTEST
    TMVAssert(A.HermOK());
#endif

    if (A.isconj()) LDL_Decompose(A.Conjugate(),xD.Conjugate(),
	P,logdet,signdet);
    else if (A.uplo() == Upper) {
      if (A.isherm()) LDL_Decompose(A.Adjoint(),xD,P,logdet,signdet);
      else LDL_Decompose(A.Transpose(),xD,P,logdet,signdet);
    } else {
      xD.Zero();
#ifdef LAP
      if (!A.iscm()) {
	Matrix<T,ColMajor> temp(A.size(),A.size());
	SymMatrixView<T> A2 = A.isherm() ?
	  HermMatrixViewOf(temp,Lower) :
	  SymMatrixViewOf(temp,Lower);
	A2 = A;
	LapLDL_Decompose(A2,xD,P,logdet,signdet);
	A = A2;
      } else {
	LapLDL_Decompose(A,xD,P,logdet,signdet);
      }
#else
      if (A.isherm()) NonLapLDL_Decompose<true>(A,xD,P,logdet,signdet);
      else NonLapLDL_Decompose<false>(A,xD,P,logdet,signdet);
#endif
      // Correct the value of logdet if det == 0
      if (signdet == T(0)) logdet = std::log(ABS(signdet));
    }

#ifdef XTEST
    TMVAssert(A.HermOK());
#endif
#ifdef XDEBUG
    //cout<<"Done LDL_Decomp\n";
    LowerTriMatrix<T,UnitDiag> L = A.LowerTri(UnitDiag);
    Matrix<T> DD(A.size(),A.size(),T(0));
    DD.diag() = A.diag();
    DD.diag(-1) = xD;
    DD.diag(1) = A.isherm() ? xD.Conjugate() : xD;
    Matrix<T> A2 = L*DD*(A.isherm() ? L.Adjoint() : L.Transpose());
    A2.ReversePermuteRows(P);
    A2.ReversePermuteCols(P);
    if (Norm(A2-A0) > 0.0001*(Norm(A0)+NormSq(L)*Norm(DD))) {
      cerr<<"LDL_Decompose\n";
      cerr<<"A0 = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"A -> "<<A<<endl;
      cerr<<"L = "<<L<<endl;
      cerr<<"D = "<<DD<<endl;
      cerr<<"A2 = "<<A2<<endl;
#ifdef LAP
      cerr<<"Compare to NonLap version:\n";
      auto_ptr<SymMatrix<T,Lower,ColMajor> > A3S(0);
      auto_ptr<HermMatrix<T,Lower,ColMajor> > A3H(0);
      auto_ptr<SymMatrixView<T> > A3(0);
      if (A.isherm()) {
	A3H.reset(new HermMatrix<T,Lower,ColMajor>(A0));
	A3.reset(new SymMatrixView<T>(A3H->View()));
      } else {
	A3S.reset(new SymMatrix<T,Lower,ColMajor>(A0));
	A3.reset(new SymMatrixView<T>(A3S->View()));
      }
      Vector<T> xD3(xD.size(),T(0));
      auto_array<int> P3(new int[A.size()]);
      RT logdet3(1);
      T signdet3(1);
      if (A.isherm()) 
	NonLapLDL_Decompose<true>(*A3,xD3.View(),P3.get(),logdet3,det3);
      else
	NonLapLDL_Decompose<false>(*A3,xD3.View(),P3.get(),logdet3,det3);
      cerr<<"A3 = "<<*A3<<endl;
      cerr<<"A = "<<A<<endl;
      cerr<<"Norm(diff) = "<<Norm(A-*A3)<<endl;
      cerr<<"xD3 = "<<xD3<<endl;
      cerr<<"xD = "<<xD<<endl;
      cerr<<"Norm(diff) = "<<Norm(xD-xD3)<<endl;
      cerr<<"P3 = ";
      for(int i=0;i<int(A.size());i++) cerr<<P3[i]<<" ";
      cerr<<"\nP  = ";
      for(int i=0;i<int(A.size());i++) cerr<<P[i]<<" ";
      cerr<<endl;
      cerr<<"det3 = "<<logdet3<<"  "<<signdet3<<endl;
      cerr<<"det = "<<logdet<<"  "<<signdet<<endl;
      abort();
#endif
    }
#endif
  }

  template <class T> void LDL_Decompose(
      const SymMatrixView<T>& A, const SymBandMatrixView<T>& D, int* P)
  {
#ifdef XDEBUG
    Matrix<T> A0 = A;
#endif
    TMVAssert(D.size() == A.size());
    TMVAssert(A.isherm() == D.isherm());
    TMVAssert(D.nlo() == 1);

    RT ld(0);
    T d(0);
    LDL_Decompose(A,D.diag(-1),P,ld,d);
    D.diag() = A.diag();
#ifdef XDEBUG
    LowerTriMatrix<T,UnitDiag> L = A.LowerTri(UnitDiag);
    Matrix<T> A2 = L * D * (A.isherm() ? L.Adjoint() : L.Transpose());
    A2.ReversePermuteRows(P);
    A2.ReversePermuteCols(P);
    //cerr<<"LDL_Decompose: A0 = "<<Type(A)<<"  "<<A0<<std::endl;
    //cerr<<"A -> "<<A<<std::endl;
    //cerr<<"D = "<<Type(D)<<"  "<<D<<std::endl;
    //cerr<<"L = "<<L<<std::endl;
    //cerr<<"LDLt = "<<L*D*(A.isherm()?L.Adjoint():L.Transpose());
    //cerr<<"D.-1 = "<<Type(D.diag(-1))<<"  "<<D.diag(-1)<<std::endl;
    //cerr<<"D.0 = "<<D.diag()<<std::endl;
    //cerr<<"D.1 = "<<D.diag(1)<<std::endl;
    //cerr<<"PLDLtPt = "<<A2<<std::endl;
    //cerr<<"cf A0 = "<<A0<<std::endl;
    if (Norm(A2-A0) > 0.001*Norm(A0)) {
      cerr<<"LDL_Decompose: A0 = "<<Type(A)<<"  "<<A0<<std::endl;
      cerr<<"A -> "<<A<<std::endl;
      cerr<<"L = "<<L<<std::endl;
      cerr<<"D = "<<D<<std::endl;
      cerr<<"LDLt = "<<L*D*(A.isherm()?L.Adjoint():L.Transpose());
      cerr<<"PLDLtPt = "<<A2<<std::endl;
      cerr<<"Norm(diff) = "<<Norm(A2-A0)<<std::endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_SymLDLDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


