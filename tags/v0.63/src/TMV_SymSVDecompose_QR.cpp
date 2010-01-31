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


#include "TMV_Blas.h"
#include "TMV_SymSVDiv.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_Matrix.h"
#include "TMV_Givens.h"

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_DiagMatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

#define RT TMV_RealType(T)

    template <class T> 
    void HermTridiagonalChopSmallElements(
        const VectorView<T>& D, const VectorView<T>& E)
    {
        // This routines sets to 0 any elements in E or D which
        // are essentially 0, given the machine precision:
        // if |E(i)| < TMV_Epsilon() * (|D(i)| + |D(i+1)|) then E(i) <- 0
        // if |D(i)| < TMV_Epsilon() * |T| then D(i) <- 0
        TMVAssert(isReal(T()));
        TMVAssert(E.size() == D.size()-1);
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        RT eps = TMV_Epsilon<T>();
        RT dthresh = eps * TMV_SQRT(NormSq(D) + 2*NormSq(E));

        T* Di = D.ptr();
        T* Ei = E.ptr();
        RT absDim1 = TMV_ABS(*Di);
        if (absDim1 < dthresh) *Di = T(0);
        ++Di;
        for(int k=E.size();k>0;--k,++Di,++Ei) {
            RT absDi = TMV_ABS(*Di);
#ifdef TMVFLDEBUG
            TMVAssert(Di >= D.first);
            TMVAssert(Di < D.last);
#endif
            if (absDi < dthresh) *Di = T(0);
            RT ethresh = eps * (absDi + absDim1);
            if (ethresh == RT(0)) ethresh = dthresh;
#ifdef TMVFLDEBUG
            TMVAssert(Ei >= E.first);
            TMVAssert(Ei < E.last);
#endif
            if (TMV_ABS(*Ei) < ethresh) *Ei = T(0);
            absDim1 = absDi;
        }
    }

    template <class T> 
    static T TridiagonalTrailingEigenValue(
        const VectorView<T>& D, const VectorView<T>& E)
    {
        // Return the Wilkinson choice for an eigenvalue of T, namely
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
        TMVAssert(isReal(T()));
        const int N = D.size();
        TMVAssert(int(E.size()) == N-1);
        TMVAssert(N > 1);

        T a = D(N-2);
        T c = D(N-1);
        T b = E(N-2);
        T bsq = b*b;
        T d = (c-a)/2;
        T dsq = d*d;
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
    static void ReduceHermTridiagonal(
        MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E)
    {
        // Reduce the superdiagonal elements of unreduced HermTridiagonal Matrix T 
        // (given by D,E) while maintaining U B Ut. 
        // Note: the input T must be unreduced - ie. all entries are non-zero.
        const int N = D.size();
        TMVAssert(N>0);
        TMVAssert(int(E.size()) == N-1);
        if (U) TMVAssert(int(U->rowsize()) == N); 
        TMVAssert(D.step()==1);
        TMVAssert(E.step()==1);

#ifdef XDEBUG
        Matrix<RT> TT(N,N,RT(0));
        TT.diag() = D;
        TT.diag(1) = TT.diag(-1) = E;
        const int M = U ? U->colsize() : 0;
        Matrix<T> A0(M,M);
        if (U) {
            A0 = *U * TT * U->adjoint();
        }
#endif

        if (N == 1) return;

        // The reduction is based on the QR algorithm to diagonalize the
        // unreduced symmetric tridiagonal matrix T.
        // The basic idea is as follows:
        // (see Golub and van Loan, chapter 8 for a derivation)
        //
        // if T is a symmetric tridiagonal matrix
        // and mu is (approximately) an eigenvalue of T
        // and the QR decomposition of T - mu I = Q R
        // Then, T' = R Q + mu I will be tridiagonal with the last 
        // subdiagonal element small.
        //
        // Note: T' = Qt (T-muI) Q + muI = Qt T Q.
        // There is a theorem then Q is unique if its first column is 
        // specified.  The first column is given by the Givens step for
        // D(0),E(0).  This step yields (for N = 6):
        //
        // G1 T = [ x x + 0 0 0 ]
        //        [ 0 x x 0 0 0 ]
        //        [ 0 x x x 0 0 ]
        //        [ 0 0 x x x 0 ]
        //        [ 0 0 0 x x x ]
        //        [ 0 0 0 0 x x ]
        //
        // G1 T G1t = [ x x + 0 0 0 ]
        //            [ x x x 0 0 0 ]
        //            [ + x x x 0 0 ]
        //            [ 0 0 x x x 0 ]
        //            [ 0 0 0 x x x ]
        //            [ 0 0 0 0 x x ]
        //
        // The + indicates the element which screws up the tri-diagonality of T'.
        // The rest of the Givens matrices will be the ones that chase this 
        // element down the matrix and away.
        //
        // Wilkinson (1968) suggested that a good choice for mu is
        // the eigenvalue of the trailing 2x2 block of T that is 
        // closer to the trailing diagonal element.
        //
        // At the end of this procedure, E(N-1) should be smaller than it was.
        // Note: This procedure works exactly if N=2.

        RT* Di = D.ptr();
        RT* Ei = E.ptr();

        RT mu = TridiagonalTrailingEigenValue(D,E);
        RT y = *Di - mu;  // = T00 - mu
        RT x = *Ei;       // = T10
        Givens<RT> G = GivensRotate(y,x);
        for(int i=0;;++i,++Di,++Ei) {
#ifdef TMVFLDEBUG
            TMVAssert(Di >= D.first);
            TMVAssert(Di < D.last);
            TMVAssert(Di+1 >= D.first);
            TMVAssert(Di+1 < D.last);
            TMVAssert(Ei >= E.first);
            TMVAssert(Ei < E.last);
#endif
            G.symMult(*Di,*(Di+1),*Ei);
            //RT Ei2 = *Ei; // = T01
            //G.mult(*Di,*Ei);
            //G.mult(Ei2,*(Di+1));
            //G.mult(*Di,Ei2);
            //G.mult(*Ei,*(++Di));
            if (U) G.mult(U->colPair(i,i+1).transpose());
            if (i==N-2) break;
            TMVAssert(x==RT(0));
#ifdef TMVFLDEBUG
            TMVAssert(Ei+1 >= E.first);
            TMVAssert(Ei+1 < E.last);
#endif
            G.mult(x,*(Ei+1));
            G = GivensRotate(*Ei,x);
        }
#ifdef XDEBUG
        if (U) {
            TT.diag() = D;
            TT.diag(-1) = TT.diag(1) = E;
            Matrix<T> A2 = *U * TT * U->adjoint();
            if (Norm(A2-A0) > 0.001*Norm(A0)) {
                cerr<<"Decompose from Tridiagonal:\n";
                cerr<<"A0 = "<<A0<<endl;
                cerr<<"A2 = "<<A2<<endl;
                cerr<<"U = "<<*U<<endl;
                cerr<<"TT = "<<TT<<endl;
                abort();
            }
        } 
#endif
    }

    template <class T> 
    void EigenFromTridiagonal_QR(
        MVP<T> U, const VectorView<RT>& D, const VectorView<RT>& E)
    {
        TMVAssert(D.size() == E.size()+1);
        if (U) {
            TMVAssert(U->rowsize() == D.size());
        }

        const int N = D.size();

        // We successively reduce the offdiagonals of T (E) to 0
        // using a sequence of Givens rotations. 
        // The reduction procedure tends to push the values up and left, so it 
        // makes sense to start at the lower right and work back up the matrix.
        // We also set to zero any very small values based on machine precision.
        // Loop invariant: all E(i) with i>=q are 0.
        // Initially q = N-1. (ie. All E(i) are potentially non-zero.)
        // When q = 0, we are done.
        HermTridiagonalChopSmallElements(D,E);
        for(int q = N-1; q>0; ) {
            if (E(q-1) == T(0)) --q;
            else {
                int p=q-1;
                while (p > 0 && (E(p-1) != T(0))) --p; 
                // Set p such that E(p-1) = 0 and all E(i) with p<=i<q are non-zero.
                if (U) {
                    ReduceHermTridiagonal<T>(
                        U->colRange(p,q+1),D.subVector(p,q+1),E.subVector(p,q));
                } else {
                    ReduceHermTridiagonal<T>(
                        0,D.subVector(p,q+1),E.subVector(p,q));
                }
                HermTridiagonalChopSmallElements(D,E);
            }
        }
    }

#undef RT

#define InstFile "TMV_SymSVDecompose_QR.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


