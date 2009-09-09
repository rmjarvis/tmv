///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2007                                                        //
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
// along with this program in the file gpl.txt.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



#include "TMV_Blas.h"
#include "TMV_Matrix.h"
#include "TMV_TriMatrix.h"
#include "TMV_SymMatrix.h"
#include "TMV_DiagMatrix.h"
#include "TMV_SymLUDiv.h"
#include "TMV_SymLUDiv_A.h"
#include "TMV_VectorArith.h"
#include "TMV_MatrixArith.h"
#include "TMV_SymMatrixArith.h"
#include <algorithm>

//#define XDEBUG

#ifdef XDEBUG
#include "TMV_TriMatrixArith.h"
#include "TMV_VIt.h"
#include <iostream>
using std::cerr;
using std::endl;
#endif

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
      std::swap(a,b);
      a/=d;
      b/=d;
      c = -c/d;
      if (dd) *dd = d;
    } else {
      T d = a*b - c*c;
      std::swap(a,b);
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
    TMVAssert(xD.step() == 1);
    const size_t N = A.size();
    static const RealType(T) alpha = (SQRT(RealType(T)(17))+1)/8;

#ifdef XTEST 
    TMVAssert(A.HermOK());
#endif
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
    T* Dj = D.ptr()+j1*D.step();
    for (size_t j=j1; j<N;) // ++j or j+=2 done below
    {
      //cerr<<"j = "<<j<<endl;
      //cerr<<"A = "<<A<<endl;
      bool seq1 = true;
      if (j == N-1) {
	// No permutation
	P[j] = j;
      } else {
	RealType(T) ajj = herm ? ABS(REAL(*Dj)) : ABS(*Dj);
	size_t p; // p is relative to j index, not absolute.
	RealType(T) apj = A.col(j,j,N).MaxAbsElement(&p);

	if (p == 0 || ajj >= alpha * apj) {
	  // No permutation
	  P[j] = j;
	} else {
	  p+=j;
	  T* Dp = D.ptr()+p*D.step();
	  RealType(T) app = herm ? ABS(REAL(*Dp)) : ABS(*Dp);
	  RealType(T) apq = A.row(p,j,p).MaxAbsElement();
	  if (p+1 < N) {
	    RealType(T) apq2 = A.col(p,p+1,N).MaxAbsElement();
	    apq = std::max(apq,apq2);
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

      //cerr<<"Before LU solving: A = "<<A<<endl;
      // Now the LU solving:
      if (seq1) {
	if (herm) {
	  RealType(T) dj = REAL(*Dj);
	  det *= dj;
	  if (dj != RealType(T)(0)) {
	    //cerr<<"0\n";
	    //cerr<<"A = "<<A<<endl;
	    A.col(j,j+1,N) /= dj;
	    //cerr<<"1\n";
	    //cerr<<"A = "<<A<<endl;
	    A.SubSymMatrix(j+1,N) -= 
	      dj*A.col(j,j+1,N)^A.col(j,j+1,N).Conjugate();
	    //cerr<<"2\n";
	    //cerr<<"A = "<<A<<endl;
	  }
	} else {
	  T dj = *Dj;
	  det *= dj;
	  if (dj != T(0)) {
	    //cerr<<"0b\n";
	    //cerr<<"A = "<<A<<endl;
	    A.col(j,j+1,N) /= dj;
	    //cerr<<"3\n";
	    //cerr<<"A = "<<A<<endl;
	    A.SubSymMatrix(j+1,N) -= dj*A.col(j,j+1,N)^A.col(j,j+1,N);
	    //cerr<<"4\n";
	    //cerr<<"A = "<<A<<endl;
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
	  if (herm) cerr<<"Herm: s==1\n";
	  else cerr<<"Sym: s==1\n";
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
	++j; Dj += D.step();
      } else {
	// Invert E:  E^-1 = [ x z* ]
	//                   [ z y  ]

	T x = *Dj;
	// Dj is at A(j,j), so
	// A(j+1,j) = *(Dj+A.stepi())
	T* Ax = Dj + A.stepi();
	T z = *(xD.ptr()+j) = *Ax;
	T y = *(Dj+=D.step());
	*Ax=T(0);
	T d;
	SymInvert_2x2<herm>(x,y,z,&d);
	if (herm) det *= REAL(d);
	else det *= d;
#ifdef XTEST
	RealType(T) NormEinv = NORM(x) + NORM(y) + 2*NORM(z);
	if (NormEinv < 1.) NormEinv = 1.;
	RealType(T) NormC = Norm(A.SubMatrix(j+2,N,j,j+2));
	if (NormC < 1.) NormC = 1.;
	RealType(T) eps = NormC*NormC*NormEinv;
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
	if (herm && IsComplex(T())) {
#ifdef XTEST
	  TMVAssert(NormInf(A.diag().Imag()) <= eps);
#endif
	  A.diag().SubVector(j+2,N).Imag().Zero();
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
	  if (herm) cerr<<"Herm: s==2\n";
	  else cerr<<"Sym: s==2\n";
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
	j+=2; Dj += D.step(); // Already did one += D.step() above
      }
    }
#ifdef XTEST 
    TMVAssert(A.HermOK());
#endif
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
    TMVAssert(xD.step()==1);
    TMVAssert(!xD.isconj());

#ifdef XDEBUG
    Matrix<T> A0(A);
#endif
#ifdef XTEST
    TMVAssert(A.HermOK());
#endif

    static const RealType(T) alpha = (SQRT(RealType(T)(17))+1)/8;
    const size_t N = A.size();

    VectorView<T> D = A.diag();
    // Only make LD if we're going to use it:
    Matrix<T,ColMajor> LD(
	N < SYM_LU_BLOCKSIZE ? 0 : N,
	N < SYM_LU_BLOCKSIZE ? 0 : SYM_LU_BLOCKSIZE+1);
    T* Dj = D.ptr();
    for (size_t j1=0; j1<N; ) {
      size_t j2 = std::min(j1+SYM_LU_BLOCKSIZE,N);

      if (j2 < N) { // On last loop we skip some steps.  See below.
	size_t j;
	T* LDjj = LD.ptr()+j1;  // = LD(j,jj)
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
	    RealType(T) ajj = ABS(*LDjj);
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

	      //RealType(T) app = ABS(LD(p,jj+1));
	      const T* LDpj = LD.cptr() + p + (j1+1)*LD.stepj();
	      RealType(T) app = ABS(*LDpj);
	      RealType(T) apq = LD.col(jj+1,j,p).MaxAbsElement();
	      if (p+1 < N) {
		RealType(T) apq2 = LD.col(jj+1,p+1,N).MaxAbsElement();
		apq = std::max(apq,apq2);
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
		*(D.ptr()+p*D.step()) = *Dj;
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
		  *(D.ptr()+p*D.step()) = *(Dj+D.step());
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
	      if (IsComplex(T())) {
		*Dj = REAL(*Dj);
	      }

	      RealType(T) dj = REAL(*Dj);
	      det *= dj;
	      if (dj != RealType(T)(0)) A.col(j,j+1,N) /= dj;
	    } else {
	      T dj = *Dj;
	      det *= dj;
	      if (dj != T(0)) A.col(j,j+1,N) /= dj;
	    }
	    ++j; 
	    Dj+=D.step();
	    LDjj += LD.stepj()+1;
	  } else {
	    // LD.cols(j,j+1) now hold L.cols(j,j+1) * E 
	    // Invert E:  E^-1 = [ x z* ]
	    //                   [ z y  ]

	    if (herm) *Dj = REAL(*LDjj);
	    else *Dj = *LDjj;
	    T x = *Dj;

	    //T z = *xDj = LD(j+1,jj);
	    ++LDjj; // Now LD(j+1,jj+1)
	    T z = *(xD.ptr()+j*xD.step()) = *LDjj;

	    Dj += D.step(); // Now D(j+1)
	    LDjj += LD.stepj(); // Now LD(j+1,jj+1)
	    if (herm) *Dj = REAL(*LDjj);
	    else *Dj = *LDjj;
	    T y = *Dj; // = D(j+1)

	    //A(j+1,j)=T(0);
	    *(A.ptr() + (j+1)*A.stepi() + j*A.stepj()) = T(0);

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

      } else NonBlockSymLU_Decompose<herm>(A,xD,P,det,j1);

#ifdef XDEBUG
      Matrix<T> L(A.LowerTri());
      L.Cols(j2,N).Zero();
      L.diag().SetAllTo(T(1));

      Matrix<T> DD(DiagMatrixViewOf(A.diag()));
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
	auto_array<size_t> P1(new size_t[A.size()]);
	T d1(1);;
	Matrix<T,ColMajor> m1(A.size(),A.size());
	SymMatrixView<T> A1 = herm ?
	  HermMatrixViewOf(m1,Lower) : SymMatrixViewOf(m1,Lower);
	A1.LowerTri() = LowerTriMatrixViewOf(A0);
	cerr<<"A1 = "<<Type(A1)<<A1<<endl;
	NonBlockSymLU_Decompose<herm>(A1.View(),xD1.View(),P1.get(),d1);
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
#ifdef XTEST
    TMVAssert(A.HermOK());
#endif
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
#ifdef INST_DOUBLE
  template <> inline void LapSymLU_Decompose(
      const SymMatrixView<double>& A, const VectorView<double>& xD,
      size_t* P, double& det)
  {
    TMVAssert(A.size()>0);
    TMVAssert(xD.size()+1 == A.size());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);
    TMVAssert(A.iscm());
    TMVAssert(!xD.isconj());

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
    const size_t Astep = A.stepj()+1;
    for(size_t i=0;i<A.size();++i,++pi,Aii+=Astep) {
      int ip = *pi LAPMINUS1;
      if (ip != int(i)) {
	if (ip >= 0) {
	  Swap(A.row(i,0,i),A.row(ip,0,i));
	  P[i] = ip;
	  det *= *Aii;
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

	  *(xD.ptr()+i*xD.step()) = b;
	  det *= a*c-b*b;
	  ++i; ++pi; // extra ++ for 2x2 case
	}
      } else {
	det *= *Aii;
	P[i] = i;
      }
    }
  }
  template <> inline void LapSymLU_Decompose(
      const SymMatrixView<std::complex<double> >& A,
      const VectorView<std::complex<double> >& xD,
      size_t* P, std::complex<double>& det)
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
    const size_t Astep = A.stepj()+1;
    for(size_t i=0;i<A.size();++i,++pi,Aii+=Astep) {
      int ip = *pi LAPMINUS1;
      if (ip != int(i)) {
	if (ip >= 0) {
	  Swap(A.row(i,0,i),A.row(ip,0,i));
	  P[i] = ip;
	  if (A.isherm()) det *= REAL(*Aii);
	  else det *= *Aii;
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
	  *(xD.ptr()+i*xD.step()) = b;
	  if (A.isherm()) {
	    double a = REAL(*Aii);
	    Aii += Astep;
	    double c = REAL(*Aii);
	    det *= a*c-NORM(b);
	  }
	  else {
	    std::complex<double> a = *Aii;
	    Aii += Astep;
	    std::complex<double> c = *Aii;
	    det *= a*c-b*b;
	  }
	  ++i; ++pi; // extra ++ for 2x2 case
	}
      } else {
	if (A.isherm()) det *= REAL(*Aii);
	else det *= *Aii;
	P[i] = i;
      }
    }
  }
#endif
#ifdef INST_FLOAT
  template <> inline void LapSymLU_Decompose(
      const SymMatrixView<float>& A, const VectorView<float>& xD,
      size_t* P, float& det)
  {
    TMVAssert(A.size()>0);
    TMVAssert(xD.size()+1 == A.size());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.uplo()==Lower);
    TMVAssert(A.iscm());
    TMVAssert(!xD.isconj());

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
    const size_t Astep = A.stepj()+1;
    for(size_t i=0;i<A.size();++i,++pi,Aii+=Astep) {
      int ip = *pi LAPMINUS1;
      if (ip != int(i)) {
	if (ip >= 0) {
	  Swap(A.row(i,0,i),A.row(ip,0,i));
	  P[i] = ip;
	  det *= *Aii;
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
	  float b = *(Aoff);
	  Aii += Astep;
	  float c = *Aii;

	  *(xD.ptr()+i*xD.step()) = b;
	  det *= a*c-b*b;
	  *Aoff = float(0);
	  ++i; ++pi; // extra ++ for 2x2 case
	}
      } else {
	det *= *Aii;
	P[i] = i;
      }
    }
  }
  template <> inline void LapSymLU_Decompose(
      const SymMatrixView<std::complex<float> >& A,
      const VectorView<std::complex<float> >& xD,
      size_t* P, std::complex<float>& det)
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
    const size_t Astep = A.stepj()+1;
    for(size_t i=0;i<A.size();++i,++pi,Aii+=Astep) {
      int ip = *pi LAPMINUS1;
      if (ip != int(i)) {
	if (ip >= 0) {
	  Swap(A.row(i,0,i),A.row(ip,0,i));
	  P[i] = ip;
	  if (A.isherm()) det *= REAL(*Aii);
	  else det *= *Aii;
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
	  *(xD.ptr()+i*xD.step()) = b;
	  if (A.isherm()) {
	    float a = REAL(*Aii);
	    Aii += Astep;
	    float c = REAL(*Aii);
	    det *= a*c-NORM(b);
	  }
	  else {
	    std::complex<float> a = *Aii;
	    Aii += Astep;
	    std::complex<float> c = *Aii;
	    det *= a*c-b*b;
	  }
	  ++i; ++pi; // extra ++ for 2x2 case
	}
      } else {
	if (A.isherm()) det *= REAL(*Aii);
	else det *= *Aii;
	P[i] = i;
      }
    }
  }
#endif 
#endif // LAP

  template <class T> void SymLU_Decompose(
      const SymMatrixView<T>& A, const VectorView<T>& xD, 
      size_t* P, T& det)
  {
    TMVAssert(A.size() > 0);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(xD.size()+1 == A.size());
#ifdef XDEBUG
    Matrix<T> A0(A);
    //cerr<<"Start SymLU_Decomp\n";
#endif
#ifdef XTEST
    TMVAssert(A.HermOK());
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

#ifdef XTEST
    TMVAssert(A.HermOK());
#endif
#ifdef XDEBUG
    //cerr<<"Done SymLU_Decomp\n";
    LowerTriMatrix<T,UnitDiag> L = A.LowerTri(UnitDiag);
    Matrix<T> DD(A.size(),A.size(),T(0));
    DD.diag() = A.diag();
    DD.diag(-1) = xD;
    DD.diag(1) = A.isherm() ? xD.Conjugate() : xD;
    Matrix<T> A2 = L*DD*(A.isherm() ? L.Adjoint() : L.Transpose());
    A2.ReversePermuteRows(P);
    A2.ReversePermuteCols(P);
    if (Norm(A2-A0) > 0.0001*(Norm(A0)+NormSq(L)*Norm(DD))) {
      cerr<<"SymLU_Decompose\n";
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
      auto_array<size_t> P3(new size_t[A.size()]);
      T det3(1);
      if (A.isherm()) 
	NonLapSymLU_Decompose<true>(*A3,xD3.View(),P3.get(),det3);
      else
	NonLapSymLU_Decompose<false>(*A3,xD3.View(),P3.get(),det3);
      cerr<<"A3 = "<<*A3<<endl;
      cerr<<"A = "<<A<<endl;
      cerr<<"Norm(diff) = "<<Norm(A-*A3)<<endl;
      cerr<<"xD3 = "<<xD3<<endl;
      cerr<<"xD = "<<xD<<endl;
      cerr<<"Norm(diff) = "<<Norm(xD-xD3)<<endl;
      cerr<<"P3 = ";
      for(size_t i=0;i<A.size();i++) cerr<<P3[i]<<" ";
      cerr<<"\nP  = ";
      for(size_t i=0;i<A.size();i++) cerr<<P[i]<<" ";
      cerr<<endl;
      cerr<<"det3 = "<<det3<<endl;
      cerr<<"det = "<<det<<endl;
      abort();
#endif
    }
#endif
  }

#define InstFile "TMV_SymLUDiv_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


