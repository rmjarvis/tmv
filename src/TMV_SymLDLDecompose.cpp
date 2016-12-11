///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2016                                                 //
// All rights reserved                                                       //
//                                                                           //
// The project is hosted at https://code.google.com/p/tmv-cpp/               //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis17 [at] gmail.                                                 //
//                                                                           //
// This software is licensed under a FreeBSD license.  The file              //
// TMV_LICENSE should have bee included with this distribution.              //
// It not, you can get a copy from https://code.google.com/p/tmv-cpp/.       //
//                                                                           //
// Essentially, you can use this software however you want provided that     //
// you include the TMV_LICENSE file in any distribution that uses it.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


//#define XDEBUG

#include "TMV_Blas.h"
#include "TMV_SymLDLDiv.h"
#include "tmv/TMV_SymLDLD.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_SymBandMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_SymMatrixArith.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_VIt.h"


#ifdef XDEBUG
#include "tmv/TMV_TriMatrixArith.h"
#include "tmv/TMV_SymBandMatrixArith.h"
#include "tmv/TMV_VIt.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

#define RT TMV_RealType(T)

#ifdef TMV_BLOCKSIZE
#define SYM_LDL_BLOCKSIZE TMV_BLOCKSIZE
#else
#define SYM_LDL_BLOCKSIZE 48
#endif

    //
    // Decompose
    //

    // TODO: This is a bit slower for row major, even with blocking, so write 
    // the row major version where A is decomposed into Lt D L rather than 
    // L D Lt.
    template <bool herm, class T> 
    static void NonBlockLDL_Decompose(
        SymMatrixView<T> A, VectorView<T> xD, 
        ptrdiff_t* P, RT& logdet, T& signdet, ptrdiff_t j1=0)
    {
        // Bunch-Kaufman algorithm for LDL Decomposition of a symmetric matrix.
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
        // Bunch and Kaufman modified the algorithm slightly to require only
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
        TMVAssert(isReal(T()) || herm == A.isherm());
        TMVAssert(herm || isComplex(T()));
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(A.ct()==NonConj);
        TMVAssert(A.uplo()==Lower);
        const ptrdiff_t N = A.size();
        const RT alpha = (TMV_SQRT(RT(17))+1)/8;

#ifdef XDEBUG
        std::cout<<"Start NonBlockLDL_Decompose\n";
        Matrix<T> A0(A);
#endif

        VectorView<T> D = A.diag();
        T* Dj = D.ptr()+j1*D.step();
        for (ptrdiff_t j=j1; j<N;) { // ++j or j+=2 done below
            //std::cout<<j<<"  ";
            bool seq1 = true;
            if (j == N-1) {
                //std::cout<<"j==N-1  ";
                // No permutation
                P[j] = j;
            } else {
                RT ajj = herm ? TMV_ABS(TMV_REAL(*Dj)) : TMV_ABS(*Dj);
                //std::cout<<ajj<<"  ";
                ptrdiff_t p; // p is relative to j index, not absolute.
                RT apj = A.col(j,j,N).maxAbsElement(&p);
                //std::cout<<apj<<"  ";
                // Check for underflow:
                if (TMV_Underflow(apj)) {
                    apj = ajj = RT(0);
                    p = 0;
                    A.col(j,j,N).setZero();
                }

                if (p == 0 || ajj >= alpha * apj) {
                    //std::cout<<"No perm (1)  ";
                    // No permutation
                    P[j] = j;
                } else {
                    p+=j;
                    T* Dp = D.ptr()+p*D.step();
                    RT app = herm ? TMV_ABS(TMV_REAL(*Dp)) : TMV_ABS(*Dp);
                    RT apq = A.row(p,j,p).maxAbsElement();
                    if (p+1 < N) {
                        RT apq2 = A.col(p,p+1,N).maxAbsElement();
                        apq = TMV_MAX(apq,apq2);
                    }
                    if (ajj*apq >= alpha * apj * apj) {
                        //std::cout<<"No perm (2)  ";
                        // No permutation
                        P[j] = j;
                    } else if (app >= alpha * apq) {
                        //std::cout<<"Perm (1)  ";
                        // Permute p diagonal into j spot
                        TMVAssert(p<A.size());
                        A.swapRowsCols(j,p);
                        P[j] = p;
                    } else {
                        //std::cout<<"Perm (2)  ";
                        // Permute pj element into j+1,j spot
                        // This also permutes pp element into j+1,j+1
                        seq1 = false;
                        P[j] = j;
                        P[j+1] = p;
                        if (p != j+1) {
                            TMVAssert(p<A.size());
                            A.swapRowsCols(j+1,p);
                        }
                    }
                }
            }

            // Now the LDL solving:
            if (seq1) {
                //std::cout<<"seq1\n";
                if (herm) {
                    RT dj = TMV_REAL(*Dj);
                    if (signdet != T(0)) {
                        if (dj == RT(0)) signdet = T(0);
                        else {
                            if (dj < RT(0)) signdet = -signdet;
                            logdet += TMV_LOG(TMV_ABS(dj));
                        }
                    }
                    if (dj != RT(0)) {
                        A.col(j,j+1,N) /= dj;
                        A.subSymMatrix(j+1,N) -= 
                            dj*A.col(j,j+1,N)^A.col(j,j+1,N).conjugate();
                    }
                } else {
                    T dj = *Dj;
                    if (signdet != T(0)) {
                        if (dj == T(0)) signdet = T(0);
                        else {
                            RT absd = TMV_ABS(dj);
                            signdet *= TMV_SIGN(dj,absd);
                            logdet += TMV_LOG(absd);
                        }
                    }
                    if (dj != T(0)) {
                        A.col(j,j+1,N) /= dj;
                        A.subSymMatrix(j+1,N) -= dj*A.col(j,j+1,N)^A.col(j,j+1,N);
                    }
                }
                ++j; Dj += D.step();
            } else {
                //std::cout<<"seq2\n";
                // Invert E:  E^-1 = [ x z* ]
                //                   [ z y  ]

                T x = *Dj;
                // Dj is at A(j,j), so
                // A(j+1,j) = *(Dj+A.stepi())
                T* Ax = Dj + A.stepi();
                T z = xD[j] = *Ax;
                T y = *(Dj+=D.step());
                //std::cout<<"x,y,z = "<<x<<','<<y<<','<<z<<std::endl;
#ifdef TMVFLDEBUG
                TMVAssert(Ax >= A._first);
                TMVAssert(Ax < A._last);
#endif
                *Ax=T(0);
                T d = SymInvert_2x2<herm>(x,y,z);
                //std::cout<<"x,y,z,d => "<<x<<','<<y<<','<<z<<','<<d<<std::endl;
                if (signdet != T(0)) {
                    if (d == T(0)) signdet = T(0);
                    else if (herm) {
                        if (TMV_REAL(d) < RT(0)) signdet = -signdet; 
                        logdet += TMV_LOG(TMV_ABS(TMV_REAL(d)));
                    } else {
                        RT absd = TMV_ABS(d);
                        signdet *= TMV_SIGN(d,absd);
                        logdet += TMV_LOG(absd);
                    }
                    //std::cout<<"logdet,signdet => "<<logdet<<','<<signdet<<std::endl;
                }
#ifdef XDEBUG
                RT normEinv = TMV_NORM(x) + TMV_NORM(y) + 2*TMV_NORM(z);
                if (normEinv < 1.) normEinv = 1.;
                RT normC = Norm(A.subMatrix(j+2,N,j,j+2));
                if (normC < 1.) normC = 1.;
                RT eps = normC*normC*normEinv;
                eps *= N*TMV_Epsilon<T>();
#endif

                if (A.iscm()) {
                    //std::cout<<"CM loop\n";
                    // Call C the current kx2 matrix which is equal to LE
                    // Need to right-multiply by E^-1 to get L
                    // A(j+2:N,j:j+2) = CE^-1
                    //
                    // Also need to update the rest of the matrix by -= L * Ct
                    // A(j+2:N,j+2:N) -= CE^-1Ct 
                    //                -= C (CE^-1)t (since E^-1 is Hermitian)
                    // A(p,q) -= C(p,0)*(L(q,0)*) + C(p,1)*(L(q,1)*)
                    //   p >= q, so loop from top down
                    const ptrdiff_t sj = A.stepj();
                    T* Cq0 = A.ptr() + (j+2) + j*sj;
                    T* Cq1 = Cq0 + sj;
                    for(ptrdiff_t q = j+2; q<N; ++q,++Cq0,++Cq1) {
                        T Lq0 = *Cq0 * x + *Cq1 * z;
                        T Lq1 = *Cq0 * (herm?TMV_CONJ(z):z) + *Cq1 * y;
                        //std::cout<<q<<"  "<<*Cq0<<','<<*Cq1<<"  "<<Lq0<<','<<Lq1<<std::endl;
                        A.col(q,q,N) -= A.col(j,q,N) * (herm?TMV_CONJ(Lq0):Lq0);
                        A.col(q,q,N) -= A.col(j+1,q,N) * (herm?TMV_CONJ(Lq1):Lq1);
#ifdef TMVFLDEBUG
                        TMVAssert(Cq0 >= A._first);
                        TMVAssert(Cq0 < A._last);
                        TMVAssert(Cq1 >= A._first);
                        TMVAssert(Cq1 < A._last);
#endif
                        *Cq0 = Lq0;
                        *Cq1 = Lq1;
                    }
                } else {
                    //std::cout<<"RM loop\n";
                    // A(j+2:N,j:j+2) = CE^-1
                    // A(j+2:N,j+2:N) -= CE^-1Ct
                    // A(p,q) -= L(p,0)*(C(q,0)*) + L(p,1)*(C(q,1)*)
                    //   p >= q, so loop from bottom up.
                    const ptrdiff_t si = A.stepi();
                    T* Cp0 = A.ptr() + (N-1)*si + j;
                    T* Cp1 = Cp0 + 1;
                    for(ptrdiff_t p = N-1; p>=j+2; --p,Cp0-=si,Cp1-=si) {
                        T Lp0 = *Cp0 * x + *Cp1 * z;
                        T Lp1 = *Cp0 * (herm?TMV_CONJ(z):z) + *Cp1 * y;
                        //std::cout<<p<<"  "<<*Cp0<<','<<*Cp1<<"  "<<Lp0<<','<<Lp1<<std::endl;
                        A.row(p,j+2,p+1) -= Lp0 * 
                            (herm ? A.col(j,j+2,p+1).conjugate() : A.col(j,j+2,p+1));
                        A.row(p,j+2,p+1) -= Lp1 *
                            (herm ? A.col(j+1,j+2,p+1).conjugate() : A.col(j+1,j+2,p+1));
#ifdef TMVFLDEBUG
                        TMVAssert(Cp0 >= A._first);
                        TMVAssert(Cp0 < A._last);
                        TMVAssert(Cp1 >= A._first);
                        TMVAssert(Cp1 < A._last);
#endif
                        *Cp0 = Lp0;
                        *Cp1 = Lp1;
                    }
                }
                if (herm && isComplex(T())) {
#ifdef XDEBUG
                    TMVAssert(NormInf(A.diag().imagPart()) <= eps);
#endif
                    A.diag().subVector(j+2,N).imagPart().setZero();
                }
                j+=2; Dj += D.step(); // Already did one += D.step() above
            }
            //std::cout<<"A -> "<<A<<std::endl;
        }
        TMVAssert(A.isHermOK());
#ifdef XDEBUG
        if (j1 == 0) {
            LowerTriMatrix<T,UnitDiag> L = A.lowerTri(UnitDiag);
            Matrix<T> DD(A.size(),A.size(),T(0));
            DD.diag() = A.diag();
            DD.diag(-1) = xD;
            DD.diag(1) = A.isherm() ? xD.conjugate() : xD;
            Matrix<T> A2 = L*DD*(A.isherm() ? L.adjoint() : L.transpose());
            A2.reversePermuteRows(P);
            A2.reversePermuteCols(P);
            std::cout<<"Done: D = "<<A.diag()<<std::endl;
            std::cout<<"     xD = "<<xD<<std::endl;
            std::cout<<"Norm(A2-A0) = "<<Norm(A2-A0)<<std::endl;
            std::cout<<"Norm(A0) = "<<Norm(A0)<<std::endl;
            std::cout<<"Norm(L) = "<<Norm(L)<<std::endl;
            std::cout<<"Norm(DD) = "<<Norm(DD)<<std::endl;
            if (!(Norm(A2-A0) < 1.e-3*(Norm(A0)+NormSq(L)*Norm(DD)))) {
                cerr<<"NonBlock LDL_Decompose\n";
                cerr<<"A0 = "<<TMV_Text(A)<<"  "<<A0<<endl;
                cerr<<"A -> "<<A<<endl;
                cerr<<"L = "<<L<<endl;
                cerr<<"D = "<<DD<<endl;
                cerr<<"A2 = "<<A2<<endl;
                abort();
            }
        }
#endif
    }

    template <bool herm, class T> 
    static void BlockLDL_Decompose(
        SymMatrixView<T> A, VectorView<T> xD, 
        ptrdiff_t* P, RT& logdet, T& signdet)
    {
        // Blocked version with level 3 calls

        TMVAssert(A.size()>0);
        TMVAssert(xD.size()+1 == A.size());
        TMVAssert(isReal(T()) || herm == A.isherm());
        TMVAssert(herm || isComplex(T()));
        TMVAssert(A.iscm() || A.isrm());
        TMVAssert(A.ct()==NonConj);
        TMVAssert(A.uplo()==Lower);

#ifdef XDEBUG
        std::cout<<"Start BlockLDL_Decompose\n";
        Matrix<T> A0(A);
#endif

        const RT alpha = (TMV_SQRT(RT(17))+1)/8;
        const ptrdiff_t N = A.size();

        VectorView<T> D = A.diag();
        Matrix<T,ColMajor> LD(N,SYM_LDL_BLOCKSIZE+1);
        T* Dj = D.ptr();
        for (ptrdiff_t j1=0; j1<N; ) {
            ptrdiff_t j2 = TMV_MIN(j1+SYM_LDL_BLOCKSIZE,N);

            if (j2 < N) { // On last loop we skip some steps.  See below.
                ptrdiff_t j;
                T* LDjj = LD.ptr()+j1;  // = LD(j,jj)
                for (j=j1; j<j2;) { // ++j or j+=2 done below
                    //std::cout<<j<<"  ";
                    bool seq1 = true;
                    ptrdiff_t jj = j-j1; // Col index in LD

                    LD.col(jj,j,N) = A.col(j,j,N);
                    LD.col(jj,j,N) -= A.subMatrix(j,N,j1,j) *
                        (herm ? LD.row(j,0,jj).conjugate() : LD.row(j,0,jj));

                    if (j == N-1) {
                        //std::cout<<"j==N-1  ";
                        // No permutation
                        P[j] = j;
                    } else {
                        RT ajj = herm ?
                            TMV_ABS(TMV_REAL(*LDjj)) :
                            TMV_ABS(*LDjj);
                        //std::cout<<ajj<<"  ";
                        ptrdiff_t p; // p is relative to j index, not absolute.
                        RT apj = LD.col(jj,j,N).maxAbsElement(&p);
                        //std::cout<<apj<<"  ";
                        // Check for underflow:
                        if (TMV_Underflow(apj)) {
                            apj = ajj = RT(0);
                            p = 0;
                            LD.col(jj,j,N).setZero();
                        }

                        if (p == 0 || ajj >= alpha * apj) {
                            //std::cout<<"No perm (1)  ";
                            // No permutation
                            P[j] = j;
                        } else {
                            p+=j;
                            //std::cout<<"p = "<<p<<"  ";

                            LD.col(jj+1,j,p) = A.col(p,j,p);
                            LD.col(jj+1,p,N) = A.col(p,p,N);
                            LD.col(jj+1,j,N) -= A.subMatrix(j,N,j1,j) * 
                                (herm ?
                                 LD.row(p,0,jj).conjugate() :
                                 LD.row(p,0,jj));

                            const T* LDpj = LD.cptr() + p + (jj+1)*LD.stepj();
                            RT app = herm ? 
                                TMV_ABS(TMV_REAL(*LDpj)) :
                                TMV_ABS(*LDpj);
                            //std::cout<<app<<"  ";
                            RT apq = LD.col(jj+1,j,p).maxAbsElement();
                            if (p+1 < N) {
                                RT apq2 = LD.col(jj+1,p+1,N).maxAbsElement();
                                apq = TMV_MAX(apq,apq2);
                            }
                            //std::cout<<apq<<"  ";
                            if (ajj*apq >= alpha * apj * apj) {
                                //std::cout<<"No perm (2)  ";
                                // No permutation
                                P[j] = j;
                            } else if (app >= alpha * apq) {
                                //std::cout<<"Perm (1)  ";
                                // Permute p diagonal into j spot
                                P[j] = p;
                                // Move updated pivot column (in LD) to j column
                                LD.col(jj,j,N) = LD.col(jj+1,j,N);
                                TMVAssert(p<A.size());
                                // Do important parts of the 
                                // A.swapRowsCols(j,p) call,
                                // We will overwrite A.col(j),
                                // so don't bother swapping into it.
                                D(p) = D(j);
                                A.row(p,j+1,p) = 
                                    herm ?
                                    A.col(j,j+1,p).conjugate() :
                                    A.col(j,j+1,p);
                                A.col(p,p+1,N) = A.col(j,p+1,N);
                                Swap(A.row(j,0,j),A.row(p,0,j));
                                // Note: this swap goes all the way back to 0,
                                // not j1
                                // Also need to do row swaps in LD
                                Swap(LD.row(j,0,jj+1),LD.row(p,0,jj+1));
                            } else {
                                //std::cout<<"Perm (2)  ";
                                // Permute pj element into j+1,j spot
                                // This also permutes pp element into j+1,j+1
                                seq1 = false;
                                P[j] = j;
                                P[j+1] = p;
                                if (p != j+1) {
                                    TMVAssert(p<A.size());
                                    // Do important parts of the 
                                    // A.swapRowsCols(j+1,p) call,
                                    // We will overwrite A.cols(j,j+1),
                                    // so don't bother 
                                    // swapping into them.
                                    D(p) = D(j+1);
                                    A.row(p,j+2,p) = 
                                        ( herm ?
                                          A.col(j+1,j+2,p).conjugate() :
                                          A.col(j+1,j+2,p) );
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
                        //std::cout<<"seq1\n";
                        // LD.col(j) now holds L.col(j)*D(j)
                        A.col(j,j,N) = LD.col(jj,j,N);
                        if (herm) {
                            if (isComplex(T())) {
#ifdef TMVFLDEBUG
                                TMVAssert(D.ptr() >= A._first);
                                TMVAssert(D.ptr() < A._last);
#endif
                                *Dj = TMV_REAL(*Dj);
                            }

                            RT dj = TMV_REAL(*Dj);
                            if (signdet != T(0)) {
                                if (dj == RT(0)) signdet = T(0);
                                else {
                                    if (dj < RT(0)) signdet = -signdet;
                                    logdet += TMV_LOG(TMV_ABS(dj));
                                }
                            }
                            if (dj != RT(0)) A.col(j,j+1,N) /= dj;
                        } else {
                            T dj = *Dj;
                            if (signdet != T(0)) {
                                if (dj == RT(0)) signdet = T(0);
                                else {
                                    RT absd = TMV_ABS(dj);
                                    signdet *= TMV_SIGN(dj,absd);
                                    logdet += TMV_LOG(absd);
                                }
                            }
                            if (dj != T(0)) A.col(j,j+1,N) /= dj;
                        }
                        ++j; 
                        Dj+=D.step();
                        LDjj += LD.stepj()+1;
                    } else {
                        //std::cout<<"seq2\n";
                        // LD.cols(j,j+1) now hold L.cols(j,j+1) * E 
                        // Invert E:  E^-1 = [ x z* ]
                        //                   [ z y  ]

#ifdef TMVFLDEBUG
                        TMVAssert(Dj >= A._first);
                        TMVAssert(Dj < A._last);
#endif
                        if (herm) *Dj = TMV_REAL(*LDjj);
                        else *Dj = *LDjj;
                        T x = *Dj;

                        ++LDjj; // Now LD(j+1,jj+1)
                        T z = xD[j] = *LDjj;

                        Dj += D.step(); // Now D(j+1)
                        LDjj += LD.stepj(); // Now LD(j+1,jj+1)
#ifdef TMVFLDEBUG
                        TMVAssert(Dj >= A._first);
                        TMVAssert(Dj < A._last);
#endif
                        if (herm) *Dj = TMV_REAL(*LDjj);
                        else *Dj = *LDjj;
                        T y = *Dj; // = D(j+1)

                        A(j+1,j)=T(0);

                        T d = SymInvert_2x2<herm>(x,y,z);
                        if (signdet != T(0)) {
                            if (d == RT(0)) signdet = T(0);
                            else if (herm) {
                                RT absd = TMV_ABS(TMV_REAL(d));
                                if (TMV_REAL(d) < 0) signdet = -signdet;
                                logdet += TMV_LOG(absd);
                            } else {
                                RT absd = TMV_ABS(d);
                                signdet *= TMV_SIGN(d,absd);
                                logdet += TMV_LOG(absd);
                            }
                        }

                        // A(j+2:N,j:j+2) = LD E^-1
                        const ptrdiff_t ldsj = LD.stepj();
                        T* LDq0 = LD.ptr() + (j+2) + jj*ldsj;
                        T* LDq1 = LDq0 + ldsj;
                        if (A.iscm()) {
                            const ptrdiff_t sj = A.stepj();
                            T* Aq0 = A.ptr() + (j+2) + j*sj;
                            T* Aq1 = Aq0 + sj;
                            for(ptrdiff_t q=j+2; q<N; ++q,++Aq0,++Aq1,++LDq0,++LDq1) {
#ifdef TMVFLDEBUG
                                TMVAssert(Aq0 >= A._first);
                                TMVAssert(Aq0 < A._last);
                                TMVAssert(Aq1 >= A._first);
                                TMVAssert(Aq1 < A._last);
#endif
                                *Aq0 = *LDq0 * x + *LDq1 * z;
                                *Aq1 = *LDq0 * (herm?TMV_CONJ(z):z) + *LDq1 * y;
                            }
                        } else {
                            const ptrdiff_t si = A.stepi();
                            T* Aq0 = A.ptr() + (j+2)*si + j;
                            T* Aq1 = Aq0 + 1;
                            for(ptrdiff_t q=j+2; q<N; ++q,Aq0+=si,Aq1+=si,++LDq0,++LDq1) {
#ifdef TMVFLDEBUG
                                TMVAssert(Aq0 >= A._first);
                                TMVAssert(Aq0 < A._last);
                                TMVAssert(Aq1 >= A._first);
                                TMVAssert(Aq1 < A._last);
#endif
                                *Aq0 = *LDq0 * x + *LDq1 * z;
                                *Aq1 = *LDq0 * (herm?TMV_CONJ(z):z) + *LDq1 * y;
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
                SymMatrixView<T> A22 = A.subSymMatrix(j2,N);
                MatrixView<T> L21 = A.subMatrix(j2,N,j1,j2);
                MatrixView<T> LD21 = LD.subMatrix(j2,N,0,j2-j1);

                // A22 -= L21 * LD21t
                SymMultMM<true>(
                    T(-1),L21,herm?LD21.adjoint():LD21.transpose(),A22);

            } else {
                NonBlockLDL_Decompose<herm>(A,xD,P,logdet,signdet,j1);
            }

            j1 = j2;
        }
        TMVAssert(A.isHermOK());
#ifdef XDEBUG
        LowerTriMatrix<T,UnitDiag> L = A.lowerTri(UnitDiag);
        Matrix<T> DD(A.size(),A.size(),T(0));
        DD.diag() = A.diag();
        DD.diag(-1) = xD;
        DD.diag(1) = A.isherm() ? xD.conjugate() : xD;
        Matrix<T> A2 = L*DD*(A.isherm() ? L.adjoint() : L.transpose());
        A2.reversePermuteRows(P);
        A2.reversePermuteCols(P);
        std::cout<<"Done: D = "<<A.diag()<<std::endl;
        std::cout<<"     xD = "<<xD<<std::endl;
        std::cout<<"Norm(A2-A0) = "<<Norm(A2-A0)<<std::endl;
        std::cout<<"Norm(A0) = "<<Norm(A0)<<std::endl;
        std::cout<<"Norm(L) = "<<Norm(L)<<std::endl;
        std::cout<<"Norm(DD) = "<<Norm(DD)<<std::endl;
        if (!(Norm(A2-A0) < 1.e-3*(Norm(A0)+NormSq(L)*Norm(DD)))) {
            cerr<<"Block LDL_Decompose\n";
            cerr<<"A0 = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"A -> "<<A<<endl;
            cerr<<"L = "<<L<<endl;
            cerr<<"D = "<<DD<<endl;
            cerr<<"A2 = "<<A2<<endl;
            cerr<<"diff = "<<Matrix<T>(A0-A2).clip(1.e-3*A0.maxAbsElement())<<endl;
            cerr<<"Norm(diff) = "<<Norm(A0-A2)<<endl;
            cerr<<"Compare to NonBlock version:\n";
            auto_ptr<SymMatrix<T,Lower|ColMajor> > A3S(0);
            auto_ptr<HermMatrix<T,Lower|ColMajor> > A3H(0);
            auto_ptr<SymMatrixView<T> > A3(0);
            if (A.isherm()) {
                A3H.reset(new HermMatrix<T,Lower|ColMajor>(A0));
                A3.reset(new SymMatrixView<T>(A3H->view()));
            } else {
                A3S.reset(new SymMatrix<T,Lower|ColMajor>(A0));
                A3.reset(new SymMatrixView<T>(A3S->view()));
            }
            Vector<T> xD3(xD.size(),T(0));
            AlignedArray<ptrdiff_t> P3(A.size());
            RT logdet3(1);
            T signdet3(1);
            if (A.isherm()) 
                NonBlockLDL_Decompose<true>(
                    *A3,xD3.view(),P3.get(),logdet3,signdet3);
            else
                NonBlockLDL_Decompose<false>(
                    *A3,xD3.view(),P3.get(),logdet3,signdet3);
            cerr<<"A3 = "<<*A3<<endl;
            cerr<<"A = "<<A<<endl;
            cerr<<"diff = "<<Matrix<T>(A-*A3).clip(1.e-3*A.maxAbsElement())<<endl;
            cerr<<"Norm(diff) = "<<Norm(A-*A3)<<endl;
            cerr<<"D3 = "<<A3->diag()<<endl;
            cerr<<"D = "<<A.diag()<<endl;
            cerr<<"Norm(diff) = "<<Norm(A.diag()-A3->diag())<<endl;
            cerr<<"xD3 = "<<xD3<<endl;
            cerr<<"xD = "<<xD<<endl;
            cerr<<"Norm(diff) = "<<Norm(xD-xD3)<<endl;
            cerr<<"P3 = ";
            for(ptrdiff_t i=0;i<A.size();i++) cerr<<P3[i]<<" ";
            cerr<<"\nP  = ";
            for(ptrdiff_t i=0;i<A.size();i++) cerr<<P[i]<<" ";
            cerr<<endl;
            cerr<<"det3 = "<<logdet3<<"  "<<signdet3<<endl;
            cerr<<"det = "<<logdet<<"  "<<signdet<<endl;
            abort();
        }
#endif
    }

    template <bool herm, class T> 
    static void NonLapLDL_Decompose(
        SymMatrixView<T> A, VectorView<T> xD, 
        ptrdiff_t* P, RT& logdet, T& signdet)
    {
        TMVAssert(A.size()>0);
        TMVAssert(xD.size()+1 == A.size());
        TMVAssert(herm == A.isherm());
        TMVAssert(A.iscm() || A.isrm());
        TMVAssert(A.ct()==NonConj);
        TMVAssert(A.uplo()==Lower);

#if 1
        if (A.size() > SYM_LDL_BLOCKSIZE) 
            BlockLDL_Decompose<herm>(A,xD,P,logdet,signdet);
        else 
#endif
            NonBlockLDL_Decompose<herm>(A,xD,P,logdet,signdet);
    }

#ifdef LAP
    template <class T> 
    static inline void LapLDL_Decompose(
        SymMatrixView<T> A, VectorView<T> xD,
        ptrdiff_t* P, RT& logdet, T& signdet)
    { 
        if (A.isherm())
            NonLapLDL_Decompose<true>(A,xD,P,logdet,signdet); 
        else
            NonLapLDL_Decompose<false>(A,xD,P,logdet,signdet); 
    }
#ifdef INST_DOUBLE
    template <> 
    void LapLDL_Decompose(
        SymMatrixView<double> A, VectorView<double> xD,
        ptrdiff_t* P, double& logdet, double& signdet)
    {
        TMVAssert(A.size()>0);
        TMVAssert(xD.size()+1 == A.size());
        TMVAssert(A.ct()==NonConj);
        TMVAssert(A.uplo()==Lower);
        TMVAssert(A.iscm());

        int n = A.size();
        int lda = A.stepj();
        AlignedArray<int> lap_p(n);
        int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = n*LAP_BLOCKSIZE;
        AlignedArray<double> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#else
        int lwork = -1;
        AlignedArray<double> work(1);
        work.get()[0] = 0.;
        LAPNAME(dsytrf) (
            LAPCM LAPCH_LO,LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1);
        lwork = int(work[0]);
        work.resize(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
        LAPNAME(dsytrf) (
            LAPCM LAPCH_LO,LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1);
#ifdef LAPNOWORK
        LAP_Results(Lap_info,"dsytrf");
#else
        LAP_Results(Lap_info,int(work[0]),n,n,lwork,"dsytrf");
#endif
        // Supposedly, sytrf is supposed to successfully complete the
        // decomposition if it encounters a singularity.  But in my 
        // experience that isn't always true.  So throw a Singular 
        // exception when it reports info > 0.
        if (Lap_info > 0) throw Singular("SymMatrix in LAPACK routine dsytrf");
            
        int* pi = lap_p.get();
        double* Aii = A.ptr();
        const int Astep = A.stepj()+1;
        for(int i=0;i<n;++i,++pi,Aii+=Astep) {
            int ip = *pi LAPMINUS1;
            if (ip != i) {
                if (ip >= 0) {
                    Swap(A.row(i,0,i),A.row(ip,0,i));
                    P[i] = ip;
                    if (signdet != 0.) {
                        if (*Aii == 0.) signdet = 0.;
                        else {
                            if (*Aii < 0.) signdet = -signdet;
                            logdet += TMV_LOG(TMV_ABS(*Aii));
                        }
                    }
                } else {
                    TMVAssert(i+1 < n);
                    TMVAssert(*(pi+1) == *pi);
                    ip = -*pi LAPMINUS1;
                    P[i] = i;
                    P[i+1] = ip;
                    if (i+1 != ip) {
                        Swap(A.row(i+1,0,i),A.row(ip,0,i));
                    }
                    double a = *Aii;
                    double* Aoff = Aii+1;
                    double b = *Aoff;
                    *Aoff = 0.;
                    Aii += Astep;
                    double c = *Aii;

                    xD[i] = b;
                    double d = a*c-b*b;
                    if (signdet != 0.) {
                        if (d == 0.) signdet = 0.;
                        else {
                            if (d < 0.) signdet = -signdet;
                            logdet += TMV_LOG(TMV_ABS(d));
                        }
                    }
                    ++i; ++pi; // extra ++ for 2x2 case
                }
            } else {
                if (signdet != 0.) {
                    if (*Aii == 0.) signdet = 0.;
                    else {
                        if (*Aii < 0.) signdet = -signdet;
                        logdet += TMV_LOG(TMV_ABS(*Aii));
                    }
                }
                P[i] = i;
            }
        }
    }
    template <> 
    void LapLDL_Decompose(
        SymMatrixView<std::complex<double> > A,
        VectorView<std::complex<double> > xD,
        ptrdiff_t* P, double& logdet, std::complex<double>& signdet)
    {
        TMVAssert(A.size()>0);
        TMVAssert(xD.size()+1 == A.size());
        TMVAssert(A.ct()==NonConj);
        TMVAssert(A.uplo()==Lower);
        TMVAssert(A.iscm());
        int n = A.size();
        int lda = A.stepj();
        AlignedArray<int> lap_p(n);
        if (A.isherm()) {
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<std::complex<double> > work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<std::complex<double> > work(1);
            work.get()[0] = 0.;
            LAPNAME(zhetrf) (
                LAPCM LAPCH_LO,LAPV(n),LAPP(A.ptr()),LAPV(lda),
                LAPP(lap_p.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1);
            lwork = int(TMV_REAL(work[0]));
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(zhetrf) (
                LAPCM LAPCH_LO,LAPV(n),LAPP(A.ptr()),LAPV(lda),
                LAPP(lap_p.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"zhetrf");
#else
            LAP_Results(Lap_info,int(TMV_REAL(work[0])),n,n,lwork,"zhetrf");
#endif
            if (Lap_info > 0) 
                throw Singular("SymMatrix in LAPACK routine zhetrf");
        } else {
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<std::complex<double> > work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<std::complex<double> > work(1);
            work.get()[0] = 0.;
            LAPNAME(zsytrf) (
                LAPCM LAPCH_LO,LAPV(n),LAPP(A.ptr()),LAPV(lda),
                LAPP(lap_p.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1);
            lwork = int(TMV_REAL(work[0]));
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(zsytrf) (
                LAPCM LAPCH_LO,LAPV(n),LAPP(A.ptr()),LAPV(lda),
                LAPP(lap_p.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"zsytrf");
#else
            LAP_Results(Lap_info,int(TMV_REAL(work[0])),n,n,lwork,"zsytrf");
#endif
            if (Lap_info > 0) 
                throw Singular("SymMatrix in LAPACK routine zsytrf");
        }

        int* pi = lap_p.get();
        std::complex<double>* Aii = A.ptr();
        const int Astep = A.stepj()+1;
        for(int i=0;i<n;++i,++pi,Aii+=Astep) {
            int ip = *pi LAPMINUS1;
            if (ip != i) {
                if (ip >= 0) {
                    Swap(A.row(i,0,i),A.row(ip,0,i));
                    P[i] = ip;
                    if (signdet != 0.) {
                        if (*Aii == 0.) signdet = 0.;
                        else if (A.isherm()) {
                            if (TMV_REAL(*Aii) < 0.) signdet = -signdet;
                            logdet += TMV_LOG(TMV_ABS(TMV_REAL(*Aii)));
                        } else {
                            double absd = TMV_ABS(*Aii);
                            signdet *= TMV_SIGN(*Aii,absd);
                            logdet += TMV_LOG(absd);
                        }
                    }
                } else {
                    TMVAssert(i+1 < n);
                    TMVAssert(*(pi+1) == *pi);
                    ip = -*pi LAPMINUS1;
                    P[i] = i;
                    P[i+1] = ip;
                    if (i+1 != ip) {
                        Swap(A.row(i+1,0,i),A.row(ip,0,i));
                    }
                    std::complex<double>* Aoff = Aii+1;
                    std::complex<double> b = *Aoff;
                    *Aoff = 0.;
                    xD[i] = b;
                    if (signdet != 0.) {
                        if (A.isherm()) {
                            double a = TMV_REAL(*Aii);
                            Aii += Astep;
                            double c = TMV_REAL(*Aii);
                            double d = a*c-TMV_NORM(b);
                            if (d == 0.) signdet = 0.;
                            else {
                                if (d < 0.) signdet = -signdet;
                                logdet += TMV_LOG(TMV_ABS(d));
                            }
                        } else {
                            std::complex<double> a = *Aii;
                            Aii += Astep;
                            std::complex<double> c = *Aii;
                            std::complex<double> d = a*c-b*b;
                            if (d == 0.) signdet = 0.;
                            else {
                                double absd = TMV_ABS(d);
                                signdet *= TMV_SIGN(d,absd);
                                logdet += TMV_LOG(absd);
                            }
                        }
                    } else { Aii += Astep; }
                    ++i; ++pi; // extra ++ for 2x2 case
                }
            } else {
                if (signdet != 0.) {
                    if (*Aii == 0.) signdet = 0.;
                    else if (A.isherm()) {
                        if (TMV_REAL(*Aii) < 0.) signdet = -signdet;
                        logdet += TMV_LOG(TMV_ABS(TMV_REAL(*Aii)));
                    } else {
                        double absd = TMV_ABS(*Aii);
                        signdet *= TMV_SIGN(*Aii,absd);
                        logdet += TMV_LOG(absd);
                    }
                }
                P[i] = i;
            }
        }
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void LapLDL_Decompose(
        SymMatrixView<float> A, VectorView<float> xD,
        ptrdiff_t* P, float& logdet, float& signdet)
    {
        TMVAssert(A.size()>0);
        TMVAssert(xD.size()+1 == A.size());
        TMVAssert(A.ct()==NonConj);
        TMVAssert(A.uplo()==Lower);
        TMVAssert(A.iscm());

        int n = A.size();
        int lda = A.stepj();
        AlignedArray<int> lap_p(n);
        int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = n*LAP_BLOCKSIZE;
        AlignedArray<float> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#else
        int lwork = -1;
        AlignedArray<float> work(1);
        work.get()[0] = 0.F;
        LAPNAME(ssytrf) (
            LAPCM LAPCH_LO,LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1);
        lwork = int(work[0]);
        work.resize(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
        LAPNAME(ssytrf) (
            LAPCM LAPCH_LO,LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1);
#ifdef LAPNOWORK
        LAP_Results(Lap_info,"ssytrf");
#else
        LAP_Results(Lap_info,int(work[0]),n,n,lwork,"ssytrf");
#endif
        if (Lap_info > 0) throw Singular("SymMatrix in LAPACK routine ssytrf");
        int* pi = lap_p.get();
        float* Aii = A.ptr();
        const int Astep = A.stepj()+1;
        for(int i=0;i<n;++i,++pi,Aii+=Astep) {
            int ip = *pi LAPMINUS1;
            if (ip != i) {
                if (ip >= 0) {
                    Swap(A.row(i,0,i),A.row(ip,0,i));
                    P[i] = ip;
                    if (signdet != 0.F) {
                        if (*Aii == 0.F) signdet = 0.F;
                        else {
                            if (*Aii < 0.F) signdet = -signdet;
                            logdet += TMV_LOG(TMV_ABS(*Aii));
                        }
                    }
                } else {
                    TMVAssert(i+1 < n);
                    TMVAssert(*(pi+1) == *pi);
                    ip = -*pi LAPMINUS1;
                    P[i] = i;
                    P[i+1] = ip;
                    if (i+1 != ip) {
                        Swap(A.row(i+1,0,i),A.row(ip,0,i));
                    }
                    float a = *Aii;
                    float* Aoff = Aii+1;
                    float b = *Aoff;
                    *Aoff = 0.F;
                    Aii += Astep;
                    float c = *Aii;

                    xD[i] = b;
                    float d = a*c-b*b;
                    if (signdet != 0.F) {
                        if (d == 0.F) signdet = 0.F;
                        else {
                            if (d < 0.F) signdet = -signdet;
                            logdet += TMV_LOG(TMV_ABS(d));
                        }
                    }
                    ++i; ++pi; // extra ++ for 2x2 case
                }
            } else {
                if (signdet != 0.F) {
                    if (*Aii == 0.F) signdet = 0.F;
                    else {
                        if (*Aii < 0.F) signdet = -signdet;
                        logdet += TMV_LOG(TMV_ABS(*Aii));
                    }
                }
                P[i] = i;
            }
        }
    }
    template <> 
    void LapLDL_Decompose(
        SymMatrixView<std::complex<float> > A,
        VectorView<std::complex<float> > xD,
        ptrdiff_t* P, float& logdet, std::complex<float>& signdet)
    {
        TMVAssert(A.size()>0);
        TMVAssert(xD.size()+1 == A.size());
        TMVAssert(A.ct()==NonConj);
        TMVAssert(A.uplo()==Lower);
        TMVAssert(A.iscm());
        int n = A.size();
        int lda = A.stepj();
        AlignedArray<int> lap_p(n);
        if (A.isherm()) {
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<std::complex<float> > work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<std::complex<float> > work(1);
            work.get()[0] = 0.F;
            LAPNAME(chetrf) (
                LAPCM LAPCH_LO,LAPV(n),LAPP(A.ptr()),LAPV(lda),
                LAPP(lap_p.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1);
            lwork = int(TMV_REAL(work[0]));
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(chetrf) (
                LAPCM LAPCH_LO,LAPV(n),LAPP(A.ptr()),LAPV(lda),
                LAPP(lap_p.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"chetrf");
#else
            LAP_Results(Lap_info,int(TMV_REAL(work[0])),n,n,lwork,"chetrf");
#endif
            if (Lap_info > 0) 
                throw Singular("SymMatrix in LAPACK routine chetrf");
        } else {
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<std::complex<float> > work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<std::complex<float> > work(1);
            work.get()[0] = 0.F;
            LAPNAME(csytrf) (
                LAPCM LAPCH_LO,LAPV(n),LAPP(A.ptr()),LAPV(lda),
                LAPP(lap_p.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1);
            lwork = int(TMV_REAL(work[0]));
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(csytrf) (
                LAPCM LAPCH_LO,LAPV(n),LAPP(A.ptr()),LAPV(lda),
                LAPP(lap_p.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"csytrf");
#else
            LAP_Results(Lap_info,int(TMV_REAL(work[0])),n,n,lwork,"csytrf");
#endif
            if (Lap_info > 0) 
                throw Singular("SymMatrix in LAPACK routine csytrf");
        }
        int* pi = lap_p.get();
        std::complex<float>* Aii = A.ptr();
        const int Astep = A.stepj()+1;
        for(int i=0;i<n;++i,++pi,Aii+=Astep) {
            int ip = *pi LAPMINUS1;
            if (ip != i) {
                if (ip >= 0) {
                    Swap(A.row(i,0,i),A.row(ip,0,i));
                    P[i] = ip;
                    if (signdet != 0.F) {
                        if (*Aii == 0.F) signdet = 0.F;
                        else if (A.isherm()) {
                            if (TMV_REAL(*Aii) < 0.F) signdet = -signdet;
                            logdet += TMV_LOG(TMV_ABS(TMV_REAL(*Aii)));
                        } else {
                            float absd = TMV_ABS(*Aii);
                            signdet *= TMV_SIGN(*Aii,absd);
                            logdet += TMV_LOG(absd);
                        }
                    }
                } else {
                    TMVAssert(i+1 < n);
                    TMVAssert(*(pi+1) == *pi);
                    ip = -*pi LAPMINUS1;
                    P[i] = i;
                    P[i+1] = ip;
                    if (i+1 != ip) {
                        Swap(A.row(i+1,0,i),A.row(ip,0,i));
                    }
                    std::complex<float>* Aoff = Aii+1;
                    std::complex<float> b = *Aoff;
                    *Aoff = 0.F;
                    xD[i] = b;
                    if (signdet != 0.F) {
                        if (A.isherm()) {
                            float a = TMV_REAL(*Aii);
                            Aii += Astep;
                            float c = TMV_REAL(*Aii);
                            float d = a*c-TMV_NORM(b);
                            if (d == 0.F) signdet = 0.F;
                            else {
                                if (d < 0.F) signdet = -signdet;
                                logdet += TMV_LOG(TMV_ABS(d));
                            }
                        } else {
                            std::complex<float> a = *Aii;
                            Aii += Astep;
                            std::complex<float> c = *Aii;
                            std::complex<float> d = a*c-b*b;
                            if (d == 0.F) signdet = 0.F;
                            else {
                                float absd = TMV_ABS(d);
                                signdet *= TMV_SIGN(d,absd);
                                logdet += TMV_LOG(absd);
                            }
                        }
                    } else { Aii += Astep; }
                    ++i; ++pi; // extra ++ for 2x2 case
                }
            } else {
                if (signdet != 0.F) {
                    if (*Aii == 0.F) signdet = 0.F;
                    else if (A.isherm()) {
                        if (TMV_REAL(*Aii) < 0.F) signdet = -signdet;
                        logdet += TMV_LOG(TMV_ABS(TMV_REAL(*Aii)));
                    } else {
                        float absd = TMV_ABS(*Aii);
                        signdet *= TMV_SIGN(*Aii,absd);
                        logdet += TMV_LOG(absd);
                    }
                }
                P[i] = i;
            }
        }
    }
#endif
#endif // LAP

    template <class T> 
    void LDL_Decompose(
        SymMatrixView<T> A, VectorView<T> xD, 
        ptrdiff_t* P, RT& logdet, T& signdet)
    {
        TMVAssert(A.size() > 0);
        TMVAssert(A.iscm() || A.isrm());
        TMVAssert(xD.size()+1 == A.size());
#ifdef XDEBUG
        Matrix<T> A0(A);
#endif
        if (A.isconj()) 
            LDL_Decompose(A.conjugate(),xD.conjugate(),P,logdet,signdet);
        else if (A.uplo() == Upper) {
            if (A.isherm()) LDL_Decompose(A.adjoint(),xD,P,logdet,signdet);
            else LDL_Decompose(A.transpose(),xD,P,logdet,signdet);
        } else {
            xD.setZero();
#ifdef LAP
            if (!A.iscm()) {
                Matrix<T,ColMajor> temp(A.size(),A.size());
                if (A.isherm() && isComplex(T())) 
                    temp.diag().imagPart().setZero();
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
            if (signdet == T(0)) logdet = TMV_LOG(TMV_ABS(signdet));
        }

        TMVAssert(A.isHermOK());
#ifdef XDEBUG
        LowerTriMatrix<T,UnitDiag> L = A.lowerTri(UnitDiag);
        Matrix<T> DD(A.size(),A.size(),T(0));
        DD.diag() = A.diag();
        DD.diag(-1) = xD;
        DD.diag(1) = A.isherm() ? xD.conjugate() : xD;
        Matrix<T> A2 = L*DD*(A.isherm() ? L.adjoint() : L.transpose());
        A2.reversePermuteRows(P);
        A2.reversePermuteCols(P);
        std::cout<<"Norm(A2-A0) = "<<Norm(A2-A0)<<std::endl;
        std::cout<<"Norm(A0) = "<<Norm(A0)<<std::endl;
        std::cout<<"Norm(L) = "<<Norm(L)<<std::endl;
        std::cout<<"Norm(DD) = "<<Norm(DD)<<std::endl;
        if (!(Norm(A2-A0) < 1.e-3*(Norm(A0)+NormSq(L)*Norm(DD)))) {
            cerr<<"LDL_Decompose\n";
            cerr<<"A0 = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"A -> "<<A<<endl;
            cerr<<"L = "<<L<<endl;
            cerr<<"D = "<<DD<<endl;
            cerr<<"A2 = "<<A2<<endl;
            cerr<<"diff = "<<Matrix<T>(A0-A2).clip(1.e-3*A0.maxAbsElement())<<endl;
            cerr<<"Norm(diff) = "<<Norm(A0-A2)<<endl;
#ifdef LAP
            cerr<<"Compare to NonLap version:\n";
            auto_ptr<SymMatrix<T,Lower|ColMajor> > A3S(0);
            auto_ptr<HermMatrix<T,Lower|ColMajor> > A3H(0);
            auto_ptr<SymMatrixView<T> > A3(0);
            if (A.isherm()) {
                A3H.reset(new HermMatrix<T,Lower|ColMajor>(A0));
                A3.reset(new SymMatrixView<T>(A3H->view()));
            } else {
                A3S.reset(new SymMatrix<T,Lower|ColMajor>(A0));
                A3.reset(new SymMatrixView<T>(A3S->view()));
            }
            Vector<T> xD3(xD.size(),T(0));
            AlignedArray<ptrdiff_t> P3(A.size());
            RT logdet3(1);
            T signdet3(1);
            if (A.isherm()) 
                NonLapLDL_Decompose<true>(
                    *A3,xD3.view(),P3.get(),logdet3,signdet3);
            else
                NonLapLDL_Decompose<false>(
                    *A3,xD3.view(),P3.get(),logdet3,signdet3);
            cerr<<"A3 = "<<*A3<<endl;
            cerr<<"A = "<<A<<endl;
            cerr<<"Norm(diff) = "<<Norm(A-*A3)<<endl;
            cerr<<"D3 = "<<A3->diag()<<endl;
            cerr<<"D = "<<A.diag()<<endl;
            cerr<<"Norm(diff) = "<<Norm(A.diag()-A3->diag())<<endl;
            cerr<<"xD3 = "<<xD3<<endl;
            cerr<<"xD = "<<xD<<endl;
            cerr<<"Norm(diff) = "<<Norm(xD-xD3)<<endl;
            cerr<<"P3 = ";
            for(ptrdiff_t i=0;i<A.size();i++) cerr<<P3[i]<<" ";
            cerr<<"\nP  = ";
            for(ptrdiff_t i=0;i<A.size();i++) cerr<<P[i]<<" ";
            cerr<<endl;
            cerr<<"det3 = "<<logdet3<<"  "<<signdet3<<endl;
            cerr<<"det = "<<logdet<<"  "<<signdet<<endl;
#endif
            abort();
        }
#endif
    }

    template <class T> 
    void LDL_Decompose(
        SymMatrixView<T> A, SymBandMatrixView<T> D, ptrdiff_t* P)
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
        LowerTriMatrix<T,UnitDiag> L = A.lowerTri(UnitDiag);
        Matrix<T> A2 = L * D * (A.isherm() ? L.adjoint() : L.transpose());
        A2.reversePermuteRows(P);
        A2.reversePermuteCols(P);
        cout<<"LDL_Decompose: A0 = "<<TMV_Text(A)<<"  "<<A0<<std::endl;
        cout<<"A -> "<<A<<std::endl;
        cout<<"D = "<<TMV_Text(D)<<"  "<<D<<std::endl;
        cout<<"L = "<<L<<std::endl;
        cout<<"LDLt = "<<L*D*(A.isherm()?L.adjoint():L.transpose());
        cout<<"D.-1 = "<<TMV_Text(D.diag(-1))<<"  "<<D.diag(-1)<<std::endl;
        cout<<"D.0 = "<<D.diag()<<std::endl;
        cout<<"D.1 = "<<D.diag(1)<<std::endl;
        cout<<"PLDLtPt = "<<A2<<std::endl;
        cout<<"cf A0 = "<<A0<<std::endl;
        if (!(Norm(A2-A0) < 0.001*Norm(A0))) {
            cerr<<"LDL_Decompose: A0 = "<<TMV_Text(A)<<"  "<<A0<<std::endl;
            cerr<<"A -> "<<A<<std::endl;
            cerr<<"L = "<<L<<std::endl;
            cerr<<"D = "<<D<<std::endl;
            cerr<<"LDLt = "<<L*D*(A.isherm()?L.adjoint():L.transpose());
            cerr<<"PLDLtPt = "<<A2<<std::endl;
            cerr<<"Norm(diff) = "<<Norm(A2-A0)<<std::endl;
            abort();
        }
#endif
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_SymLDLDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


