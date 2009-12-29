///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2009                                                 //
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


//#define XDEBUG


#include "TMV_SVDiv.h"
#include "tmv/TMV_SVD.h"
#include "TMV_Givens.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_Matrix.h"

#ifdef XDEBUG
#include <iostream>
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_DiagMatrixArith.h"
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

#define RT TMV_RealType(T)

    //
    // DecomposeFromBidiagonal: QR method
    //

    template <class T> 
    void BidiagonalChopSmallElements(
        const VectorView<T>& D, const VectorView<T>& E, bool* zd)
    {
        // This routines sets to 0 any elements in D,E which
        // are essentially 0, given the machine precision:
        // if D(i)*TMV_Epsilon() == 0, the D(i) = 0
        // if E(i)*TMV_Epsilon() < abs(D(i)), the D(i) = 0
        TMVAssert(E.size() == D.size()-1);
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        RT eps = TMV_Epsilon<T>();

        T* Di = D.ptr();
        T* Ei = E.ptr();
#ifdef TMVFLDEBUG
        TMVAssert(Di >= D.first);
        TMVAssert(Di < D.last);
#endif
        if (*Di*eps == T(0)) { if(zd) *zd = true; }
        ++Di;
        for(int k=E.size();k>0;--k,++Di,++Ei) {
#ifdef TMVFLDEBUG
            TMVAssert(Di >= D.first);
            TMVAssert(Di < D.last);
#endif
            if (*Di*eps == RT(0)) { if(zd) *zd = true; }
#ifdef TMVFLDEBUG
            TMVAssert(Ei >= E.first);
            TMVAssert(Ei < E.last);
#endif
            if (TMV_MAXABS(*Ei) <= eps * TMV_MIN(TMV_MAXABS(*Di),TMV_MAXABS(*(Di-1)))) *Ei = T(0);
        }
    }

    template <class T> 
    void BidiagonalZeroFirstRow(
        MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E)
    {
        // Input D,E form an N+1 x N bidiagonal matrix
        // (eg. for N = 4)
        //     [ x 0 0 0 ]
        //     [ x x 0 0 ]
        // B = [ 0 x x 0 ]
        //     [ 0 0 x x ]
        //     [ 0 0 0 x ]
        // Zero out the first row maintaining the constancy of U B
        // using Givens transformations.
        TMVAssert(E.size() == D.size());
        if (U) TMVAssert(U->rowsize() == D.size()+1);
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);

        const int N = D.size();
        RT* Di = D.ptr();
        RT* Ei = E.ptr();

        RT x = *Ei;
        if (x != RT(0)) {
#ifdef TMVFLDEBUG
            TMVAssert(Ei >= E.first);
            TMVAssert(Ei < E.last);
#endif
            *Ei = RT(0);
            ++Ei;
            // Loop Invariant: x = B(0,i)
            for(int i=0; i<N; ++i,++Di,++Ei) {
#ifdef TMVFLDEBUG
                TMVAssert(Di >= D.first);
                TMVAssert(Di < D.last);
#endif
                Givens<RT> G = GivensRotate(*Di,x);
                // Make new B = G B
                if (i<N) {
#ifdef TMVFLDEBUG
                    TMVAssert(Ei >= E.first);
                    TMVAssert(Ei < E.last);
#endif
                    G.mult(*Ei,x);
                }
                // Make new U = U Gt
                if (U) G.conjMult(U->colPair(i+1,0).transpose());
            }
        }
    }

    template <class T> 
    void BidiagonalZeroLastCol(
        const VectorView<RT>& D, const VectorView<RT>& E, MVP<T> V)
    {
        // Input D,E form an N x N+1 bidiagonal matrix
        // (eg. for N = 4)
        //     [ x x 0 0 0 ]
        // B = [ 0 x x 0 0 ]
        //     [ 0 0 x x 0 ]
        //     [ 0 0 0 x x ]
        // Zero out the last col maintaining the constancy of B V
        // using Givens transformations.
        TMVAssert(E.size() == D.size());
        if (V) TMVAssert(V->colsize() == D.size()+1);
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);

        const int N = D.size();
        RT* Di = D.ptr()+N-1;
        RT* Ei = E.ptr()+N-1;

        RT x = *Ei;
        if (x != RT(0)) {
#ifdef TMVFLDEBUG
            TMVAssert(Ei >= E.first);
            TMVAssert(Ei < E.last);
#endif
            *Ei = RT(0);
            // Loop Invariant: x = B(i,N-1)
            for(int i=N-1; i>=0; --i,--Di) {
#ifdef TMVFLDEBUG
                TMVAssert(Di >= D.first);
                TMVAssert(Di < D.last);
#endif
                Givens<RT> G = GivensRotate(*Di,x);
                // Make new B = B GT
                if (i>0) {
#ifdef TMVFLDEBUG
                    TMVAssert(Ei-1 >= E.first);
                    TMVAssert(Ei-1 < E.last);
#endif
                    G.mult(*(--Ei),x);
                }
                // Make new V = G* V 
                if (V) G.conjMult(V->rowPair(i,N));
            }
        }
    }

    template <class T> 
    static RT BidiagonalTrailingEigenValue(
        const VectorView<T>& D, const VectorView<T>& E)
    {
        // Return the Wilkinson choice for an eigenvalue of T = BtB, namely
        // the eigenvalue of the trailing 2x2 block of T which is closer
        // to the last diagonal element of T.
        //
        // Trailing 2x2 block =  (Use i = N-2, j = N-1)
        // [ a  b ] = [ |Di|^2 + |Ei-1|^2      Di Ei      ]
        // [ b* c ]   [     Di* Ei*       |Dj|^2 + |Ei|^2 ]
        // 
        // mu = c - d +- sqrt(d^2+|b|^2), where d = (c-a)/2
        // if d>0 we use +, if d<0 we use -.
        // 
        // For stability when |b| is small, we rearrange this to:
        // mu = c + |b|^2/(d +- sqrt(d^2+|b|^2))
        //    = c +- |b|^2/|d|(1 + sqrt(1+|b|^2/d^2))
        const int N = D.size();
        TMVAssert(int(E.size()) == N-1);
        TMVAssert(N > 1);

        RT a = TMV_NORM(D(N-2)) + (N>2 ? TMV_NORM(E(N-3)) : RT(0));
        RT c = TMV_NORM(D(N-1)) + TMV_NORM(E(N-2));
        T b = D(N-2)*E(N-2);
        RT bsq = TMV_NORM(b);
        RT d = (c-a)/2;
        RT dsq = TMV_SQR(d);
        if (dsq < bsq) {
            RT x = TMV_SQRT(dsq+bsq);
            if (d > 0) return c - d + x;
            else return c - d - x;
        } else {
            RT x = bsq/TMV_ABS(d)/(1 + TMV_SQRT(1+bsq/dsq));
            if (d > 0) return c + x;
            else return c - x;
        }
    }

    template <class T> 
    static void ReduceBidiagonal22(
        MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E, MVP<T> V)
    {
        // Find Givens rotations which diagonalize 2x2 matrix exactly:
        //
        // [  c1 s1 ] [ d0 e ] [  c2 s2 ] = [  c1 d0    c1 e + s1 d1 ] [  c2 s2 ]
        // [ -s1 c1 ] [ 0 d1 ] [ -s2 c2 ]   [ -s1 d0   -s1 e + c1 d1 ] [ -s2 c2 ]
        // = [ c1 c2 d0 - c1 s2 e - s1 s2 d1     c1 s2 d0 + c1 c2 e + s1 c2 d1 ]
        //   [ -s1 c2 d0 + s1 s2 e - c1 s2 d1   -s1 s2 d0 - s1 c2 e + c1 c2 d1 ]
        //
        // lower left = 0 => s1 (s2 e - c2 d0) = c1 s2 d1
        //                   t1 = t2 d1 / (t2 e - d0)
        //
        // upper right = 0 => s1 c2 d1 = -c1 s2 d0 - c1 c2 e
        //                    t1 = -(t2 d0 + e) / d1
        // 
        // Equate these to get:
        // t2 d1^2 = -(t2 d0 + e)(t2 e - d0) = -t2^2 d0 e - t2 e^2 + d0^2 t2 + e d0
        // d0 e t2^2 + (d1^2 + e^2 - d0^2) t2 - e d0 = 0
        // t2^2 + ((c-a)/b) t2 - 1 = 0
        // where a = d0^2, b = d0 e, c = d1^2 + e^2 
        // are the elements of B^T B.  (c.f. BidiagonalTrailingEigenvalue)
        //
        // The solution is 
        //    t2 = b / ( d +- sqrt(d^2 + b^2) )
        //       = (-d +- sqrt(d^2 + b^2)) / b
        // where d = (c-a)/2 and we choose the smaller solution.
        //
        // From this, c2 = 1 / sqrt(1+t2^2)  and  s2 = t2 c2
        // Then t1 = -(t2 d0 + e) / d1
        // c1 = 1 / sqrt(1+t1^2) and s1 = t1 c1
        // 
        // D(0) <- d0 c1/c2
        // D(1) <- d1 c2/c1
        //
        // However, these formulae are not stable for the typically small
        // values of b and (sometimes) d.
        //
        // We use the following formulae which are more stable.
        // The derivations are messy, so we leave them as "an exercise 
        // for the interested reader". :)
        //
        // Let x = sqrt(d^2+b^2) if |d1| >= |d0|
        //     or -sqrt(d^2+b^2) if |d1| <  |d0|
        //
        // s2 = b / (x sqrt(2(1+d/x)))
        // s1 = -s2 sign(d0 d1) sqrt(1+(x+d)/a); // if x>0
        // s1 = -s2 d1/d0 sqrt(1-xmd/a);         // if x<0
        // D(0) -= d0 xmd / (a + sqrt(a(a-xmd)))
        // D(1) += d1 xmd / (a-xmd + sqrt(a(a-xmd)))
        //         where xmd = x-d = b^2/(x+d)
        //

        TMVAssert(D.size() == 2);
        TMVAssert(E.size() == 1);
        if (U) TMVAssert(U->rowsize() == 2);
        if (V) TMVAssert(V->colsize() == 2);

        RT d0 = D(0);
        RT d1 = D(1);
        RT e = E(0);

        RT a = d0*d0;
        RT b = d0*e;
        RT d = ((d1-d0)*(d1+d0)+e*e)/RT(2);
        RT x = TMV_SQRT(d*d+b*b);
        if (TMV_ABS(d1) < TMV_ABS(d0)) x = -x;
        RT xmd = b*b/(x+d);

        RT s2 = b / (x * TMV_SQRT(RT(2)*(RT(1)+d/x)));
        RT s1;
        if (x > RT(0)) {
            s1 = -s2 * TMV_SQRT(RT(1)+(x+d)/a);
            if (d0*d1 < 0) s1 = -s1;
        } else {
            s1 = -s2 * d1/(d0*TMV_SQRT(RT(1)-xmd/a));
        }

        RT c1 = TMV_SQRT(RT(1)-s1*s1);
        RT c2 = TMV_SQRT(RT(1)-s2*s2);

        // For small s, improve calculation of c:
        // s^2 = 1-c^2 = (1-c)(1+c)
        // c = 1 - s^2/(1+c)
        if (TMV_ABS(s1) < RT(1.e-2)) c1 = RT(1) - s1*s1/(RT(1)+c1);
        if (TMV_ABS(s2) < RT(1.e-2)) c2 = RT(1) - s2*s2/(RT(1)+c2);

#ifdef XDEBUG
        Matrix<RT> B(2,2); B(0,0) = d0; B(0,1) = e; B(1,0) = RT(0); B(1,1) = d1;
        Matrix<RT> g1(2,2); g1(0,0) = c1; g1(0,1) = s1; g1(1,0) = -s1; g1(1,1) = c1;
        Matrix<RT> g2(2,2); g2(0,0) = c2; g2(0,1) = s2; g2(1,0) = -s2; g2(1,1) = c2;
        Matrix<RT> S = g1 * B * g2;
        Matrix<T> A(U&&V ? U->colsize() : 0, U&&V ? V->rowsize() : 0);
        if (U && V) A = *U * B * *V;
#endif

        RT temp = a+TMV_SQRT(a*(a-xmd));
        D(0) -= d0*xmd/temp;
        D(1) += d1*xmd/(temp-xmd);
        E(0) = RT(0);

        if (U) {
            Givens<RT> G1(c1,s1);
            G1.mult(U->transpose());
        }
        if (V) {
            Givens<RT> G2(c2,-s2);
            G2.mult(*V);
        }
#ifdef XDEBUG
        if (U && V) {
            B.diag() = D;
            B.diag(1) = E;
            Matrix<T> A2 = *U * B * *V;
            //cout<<"B = "<<B<<endl;
            //cout<<"V = "<<*V<<endl;
            //cout<<"UBV = "<<A2<<endl;
            //cout<<"Done 22: Norm(A2-A) = "<<Norm(A2-A)<<endl;
            if (Norm(A2-A) > 0.001*Norm(A)) {
                cerr<<"ReduceBidiagonal22\n";
                cerr<<"B = "<<B<<endl;
                cerr<<"A = "<<A<<endl;
                cerr<<"S = "<<S<<endl;
                cerr<<"G1 = "<<g1<<endl;
                cerr<<"G2 = "<<g2<<endl;
                cerr<<"A2 = "<<A2<<endl;
                cerr<<"Norm(diff) = "<<Norm(A2-A)<<endl;
                abort();
            }
        }
#endif
    }

    template <class T> 
    static void ReduceBidiagonal(
        MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E, MVP<T> V)
    {
#ifdef XDEBUG
        //cout<<"Start Reduce Bidiagonal QR:\n";
        //cout<<"D = "<<D<<endl;
        //cout<<"E = "<<E<<endl;
        Matrix<RT> B(D.size(),D.size(),RT(0));
        Vector<RT> D0 = D;
        Vector<RT> E0 = E;
        B.diag() = D0;
        B.diag(1) = E0;
        Matrix<T> A0(U&&V ? U->colsize() : D.size(),
                     U&&V ? V->rowsize() : D.size());
        if (U && V) A0 = (*U) * B * (*V);
        else A0 = B;
        //cout<<"A0 = "<<A0<<endl;
#endif
        // Reduce the superdiagonal elements of Bidiagonal Matrix B 
        // (given by D,E) while maintaining U B V. 
        // Note: the input B must be unreduced - ie. all entries are non-zero.
        const int N = D.size();
        TMVAssert(N>0);
        TMVAssert(int(E.size()) == N-1);
        if (U) TMVAssert(int(U->rowsize()) == N);
        if (V) TMVAssert(int(V->colsize()) == N);
        TMVAssert(D.step()==1);
        TMVAssert(E.step()==1);
        TMVAssert(D.minAbsElement() > RT(0));
        TMVAssert(E.minAbsElement() > RT(0));

        if (N == 1) return;
        if (N == 2) return ReduceBidiagonal22(U,D,E,V);

        // The reduction is based on the QR algorithm to diagonalize the
        // unreduced symmetric tridiagonal matrix T = BtB
        // The basic idea is as follows:
        // (see Golub and van Loan, chapter 8 for a derivation)
        //
        // if T is a symmetric tridiagonal matrix
        // and mu is (approximately) an eigenvalue of T
        // and the QR decomposition of T - mu I = V R
        // Then, T' = R V + mu I will be tridiagonal with the last 
        // subdiagonal element small.
        // (Note: T' = Vt (T-muI) V + muI = Vt T V.)
        //
        // Wilkinson (1968) suggested that a good choice for mu is
        // the eigenvalue of the trailing 2x2 block of T that is 
        // closer to the trailing diagonal element.
        //
        // Rather than explicitly forming T = BtB and doing this
        // procedure, Golub and van Load show that it can be done
        // in place.
        // If T' = Vt T V, 
        // then Bt'B' = Vt Bt B V
        // B' = U B V for some U
        //
        // So, start with the first Givens matrix in the QR algorithm for T:
        // G0 [ T00 - mu ] = [ x ]
        //    [   T10    ]   [ 0 ]
        // We apply this to the right of B which has the effect: (for N=5)
        //              [ x x 0 0 0 ]
        //              [ + x x 0 0 ]
        // B <- B G0T = [ 0 0 x x 0 ]
        //              [ 0 0 0 x x ]
        //              [ 0 0 0 0 x ]
        // The + is the element which screws up the bidiagonal structure.
        // The rest of the procedure simply involves chasing this + down
        // the diagonal using Givens rotations.
        // For each Givens rotation we use, we also Multiply U or V by the
        // adjoint to maintain the constancy of U B V.
        //
        // At the end of this procedure, E(N-1) should be smaller than it was.
        // Note: This procedure works exactly if N=2.
        RT* Di = D.ptr();
        RT* Ei = E.ptr();

        RT mu = BidiagonalTrailingEigenValue(D,E);
        RT y = TMV_NORM(*Di) - mu;  // = T00 - mu
        RT x = TMV_CONJ(*Di)*(*Ei);  // = T10
        Givens<RT> G = GivensRotate(y,x);
        for(int i=1;i<N;++i) {
            G.mult(*Di,*Ei);
            if (V) G.conjMult(V->rowPair(i-1,i));
            TMVAssert(x==RT(0));
            G.mult(x,*(++Di)); // x = B(i,i-1)
            G = GivensRotate(*(Di-1),x);
            G.mult(*Ei,*Di);
            if (U) G.conjMult(U->colPair(i-1,i).transpose());
            if (i < N-1) {
                TMVAssert(x==RT(0));
                G.mult(x,*(++Ei)); // x = B(i-1,i+1)
                G = GivensRotate(*(Ei-1),x);
            }
        }
#ifdef XDEBUG
        if (U && V) {
            B.diag() = D;
            B.diag(1) = E;
            Matrix<T> AA = (*U) * B * (*V);
            if (Norm(A0-AA) > 0.001*Norm(A0)) {
                cerr<<"ReduceBidiagonal: \n";
                cerr<<"input D = "<<D0<<endl;
                cerr<<"input E = "<<E0<<endl;
                cerr<<"UBV = "<<A0<<endl;
                cerr<<"UBV => "<<AA<<endl;
                cerr<<"U = "<<*U<<endl;
                cerr<<"D = "<<D<<endl;
                cerr<<"E = "<<E<<endl;
                cerr<<"V = "<<*V<<endl;
                abort();
            }
        }
        if (Norm(D-D0) + Norm(E-E0) == RT(0)) {
            cerr<<"ReduceBidiagonal didn't reduce:\n";
            cerr<<"input D = "<<D0<<endl;
            cerr<<"input E = "<<E0<<endl;
            cerr<<"output D = "<<D<<endl;
            cerr<<"output E = "<<E<<endl;
            abort();
        }
#endif
    }

    template <class T> 
    void SV_DecomposeFromBidiagonal_QR(
        MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E, MVP<T> V)
    {
#ifdef XDEBUG
        //cout<<"Start Decompose from Bidiagonal QR:\n";
        //if (U) cout<<"U = "<<*U<<endl;
        //if (V) cout<<"V = "<<*V<<endl;
        //cout<<"D = "<<D<<endl;
        //cout<<"E = "<<E<<endl;
        Matrix<RT> B(D.size(),D.size(),RT(0));
        B.diag() = D;
        B.diag(1) = E;
        Matrix<T> A0(U&&V ? U->colsize() : D.size(),
                     U&&V ? V->rowsize() : D.size());
        if (U && V) A0 = (*U) * B * (*V);
        else A0 = B;
        //cout<<"A0 = "<<A0<<endl;
#endif


        if (E.size() == 0) return;

        TMVAssert(D.size() == E.size()+1);
        if (U) TMVAssert(U->rowsize() == D.size());
        if (V) TMVAssert(V->colsize() == D.size());
        TMVAssert(D.minAbsElement() > RT(0));
        TMVAssert(E.minAbsElement() > RT(0));

        // We successively reduce the superdiagonal of B (E) to 0
        // using a sequence of Givens rotations. 
        // The reduction procedure tends to push the values up and left, so it 
        // makes sense to start at the lower right and work back up the matrix.
        // We also set to zero any very small values based on machine precision.
        // Loop invariant: all E(i) with i>=q are 0.
        // Initially q = N-1. (ie. All E(i) are potentially non-zero.)
        // When q = 0, we are done.

        for(int q = E.size(); q>0; ) {
            if (E(q-1) == T(0)) --q;
            else {
                int p=q-1;
                while (p > 0 && (E(p-1) != T(0))) --p;
                // Set p such that E(p-1) = 0 and all E(i) with p<=i<q are non-zero.

                if (U) {
                    if (V) ReduceBidiagonal<T>(
                        U->colRange(p,q+1),D.subVector(p,q+1),
                        E.subVector(p,q),V->rowRange(p,q+1));
                    else ReduceBidiagonal<T>(
                        U->colRange(p,q+1),D.subVector(p,q+1),
                        E.subVector(p,q),0);
                } else {
                    if (V) ReduceBidiagonal<T>(
                        0,D.subVector(p,q+1),
                        E.subVector(p,q),V->rowRange(p,q+1));
                    else ReduceBidiagonal<T>(
                        0,D.subVector(p,q+1),E.subVector(p,q),0);
                }

                bool newzeroD = false;
                BidiagonalChopSmallElements(
                    D.subVector(p,q+1),E.subVector(p,q),&newzeroD);
                // Check that we haven't introduced new 0's in the D vector.
                // If we have, we need to go back to the original calling 
                // routine which deals with them.
                if (newzeroD) {
                    if (U) {
                        if (V) DoSVDecomposeFromBidiagonal<T>(
                            U->colRange(p,q+1),D.subVector(p,q+1),
                            E.subVector(p,q),V->rowRange(p,q+1),false,false);
                        else DoSVDecomposeFromBidiagonal<T>(
                            U->colRange(p,q+1),D.subVector(p,q+1),
                            E.subVector(p,q),0,false,false);
                    } else {
                        if (V) DoSVDecomposeFromBidiagonal<T>(
                            0,D.subVector(p,q+1),
                            E.subVector(p,q),V->rowRange(p,q+1),false,false);
                        else DoSVDecomposeFromBidiagonal<T>(
                            0,D.subVector(p,q+1),
                            E.subVector(p,q),0,false,false);
                    }
                    q = p;
                }
            }
        }
#ifdef XDEBUG
        if (U && V) {
            Matrix<T> AA = (*U) * DiagMatrixViewOf(D) * (*V);
            cout<<"Done QR Norm(A0-AA) = "<<Norm(A0-AA)<<endl;
            cout<<"Norm(UtU-1) = "<<Norm(U->adjoint()*(*U)-T(1))<<endl;
            cout<<"Norm(VtV-1) = "<<Norm(V->adjoint()*(*V)-T(1))<<endl;
            cout<<"Norm(VVt-1) = "<<Norm((*V)*V->adjoint()-T(1))<<endl;
            if (Norm(A0-AA) > 0.001*Norm(A0)) {
                cerr<<"SV_DecomposeFromBidiagonal QR: \n";
                cerr<<"input B = "<<B<<endl;
                cerr<<"UBV = "<<A0<<endl;
                cerr<<"USV = "<<AA<<endl;
                cerr<<"U = "<<*U<<endl;
                cerr<<"S = "<<D<<endl;
                cerr<<"V = "<<*V<<endl;
                abort();
            }
        }
#endif
    }

#undef RT

#define InstFile "TMV_SVDecompose_QR.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


