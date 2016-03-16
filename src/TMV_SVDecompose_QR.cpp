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


#include "TMV_SVDiv.h"
#include "tmv/TMV_SVD.h"
#include "tmv/TMV_Givens.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_Matrix.h"
#include <iostream>

#ifdef XDEBUG
#define THRESH 1.e-15
#include <iostream>
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_DiagMatrixArith.h"
#include "tmv/TMV_VIt.h"
#define dbgcout std::cout
//#define dbgcout if(false) std::cout
using std::cerr;
using std::endl;
#else
#define dbgcout if(false) std::cout
#endif

namespace tmv {

#define RT TMV_RealType(T)

    //
    // DecomposeFromBidiagonal: QR method
    //

    template <class T> 
    void BidiagonalChopSmallElements(
        VectorView<T> D, VectorView<T> E, bool* zd)
    {
        // This routines sets to 0 any elements in D,E which
        // are essentially 0, given the machine precision:
        // if |D(i)|^2*Epsilon == 0, the D(i) = 0
        // if |E(i)| < Epsilon * (|D(i)| + |D(i+1)|), then E(i) = 0
        TMVAssert(E.size() == D.size()-1);
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        const RT eps = TMV_Epsilon<T>();

        T* Di = D.ptr();
        T* Ei = E.ptr();
#ifdef TMVFLDEBUG
        TMVAssert(Di >= D._first);
        TMVAssert(Di < D._last);
#endif
        if (TMV_Underflow(TMV_NORM(*Di))) {
            *Di = T(0);
            if(zd) *zd = true; 
        }
        ++Di;
        for(ptrdiff_t k=E.size();k>0;--k,++Di,++Ei) {
#ifdef TMVFLDEBUG
            TMVAssert(Di >= D._first);
            TMVAssert(Di < D._last);
#endif
            if (TMV_Underflow(TMV_NORM(*Di))) {
                *Di = T(0);
                if(zd) *zd = true; 
            }
#ifdef TMVFLDEBUG
            TMVAssert(Ei >= E._first);
            TMVAssert(Ei < E._last);
#endif
            // Do it as !(|e| > eps (|d1| + |d2|)) so nan's will get set to 
            // zero too and not iterate forever.
            if ( !(TMV_MAXABS(*Ei) > eps*(TMV_MAXABS(*Di)+TMV_MAXABS(*(Di-1))))
                 || TMV_Underflow(*Ei) ) {
                *Ei = T(0);
            }

            // This last one checks that when the product of a D*E underflows
            // that at least one of the two values have been set to zero.
            // Otherwise, set the smaller one to 0.
            // e.g. if underflow happens at 1.e-310, and you have
            // Di = 1.e-154, Ei = 1.e-157, then none of the above tests
            // will trigger a zero.  But the product will underflow, leading
            // to problems with Reduce not reducing.
            if (TMV_Underflow(*Di * *Ei) && *Di!=T(0) && *Ei!=T(0)) {
                if (TMV_MAXABS(*Ei) <= TMV_MAXABS(*Di) ) *Ei = T(0);
                else *Di = T(0);
            }
            if (TMV_Underflow(*(Di-1) * *Ei) && *(Di-1)!=T(0) && *Ei!=T(0)) {
                if (TMV_MAXABS(*Ei) <= TMV_MAXABS(*(Di-1)) ) *Ei = T(0);
                else *(Di-1) = T(0);
            }
        }
    }

    template <class T> 
    void BidiagonalZeroFirstRow(
        MatrixView<T> U, VectorView<RT> D, VectorView<RT> E)
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
        if (U.cptr()) TMVAssert(U.rowsize() == D.size()+1);
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);

        const ptrdiff_t N = D.size();
        RT* Di = D.ptr();
        RT* Ei = E.ptr();

        RT x = *Ei;
        if (x != RT(0)) {
#ifdef TMVFLDEBUG
            TMVAssert(Ei >= E._first);
            TMVAssert(Ei < E._last);
#endif
            *Ei = RT(0);
            ++Ei;
            // Loop Invariant: x = B(0,i)
            for(ptrdiff_t i=0; i<N; ++i,++Di,++Ei) {
#ifdef TMVFLDEBUG
                TMVAssert(Di >= D._first);
                TMVAssert(Di < D._last);
#endif
                Givens<RT> G = GivensRotate(*Di,x);
                // Make new B = G B
                if (i<N) {
#ifdef TMVFLDEBUG
                    TMVAssert(Ei >= E._first);
                    TMVAssert(Ei < E._last);
#endif
                    G.mult(*Ei,x);
                }
                // Make new U = U Gt
                if (U.cptr()) G.conjMult(U.colPair(i+1,0).transpose());
            }
        }
    }

    template <class T> 
    void BidiagonalZeroLastCol(
        VectorView<RT> D, VectorView<RT> E, MatrixView<T> Vt)
    {
        // Input D,E form an N x N+1 bidiagonal matrix
        // (eg. for N = 4)
        //     [ x x 0 0 0 ]
        // B = [ 0 x x 0 0 ]
        //     [ 0 0 x x 0 ]
        //     [ 0 0 0 x x ]
        // Zero out the last col maintaining the constancy of B Vt
        // using Givens transformations.
        TMVAssert(E.size() == D.size());
        if (Vt.cptr()) TMVAssert(Vt.colsize() == D.size()+1);
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);

        const ptrdiff_t N = D.size();
        RT* Di = D.ptr()+N-1;
        RT* Ei = E.ptr()+N-1;

        RT x = *Ei;
        if (x != RT(0)) {
#ifdef TMVFLDEBUG
            TMVAssert(Ei >= E._first);
            TMVAssert(Ei < E._last);
#endif
            *Ei = RT(0);
            // Loop Invariant: x = B(i,N-1)
            for(ptrdiff_t i=N-1; i>=0; --i,--Di) {
#ifdef TMVFLDEBUG
                TMVAssert(Di >= D._first);
                TMVAssert(Di < D._last);
#endif
                Givens<RT> G = GivensRotate(*Di,x);
                // Make new B = B GT
                if (i>0) {
#ifdef TMVFLDEBUG
                    TMVAssert(Ei-1 >= E._first);
                    TMVAssert(Ei-1 < E._last);
#endif
                    G.mult(*(--Ei),x);
                }
                // Make new Vt = G* Vt 
                if (Vt.cptr()) G.conjMult(Vt.rowPair(i,N));
            }
        }
    }

    template <class T> 
    static RT BidiagonalTrailingEigenValue(
        VectorView<T> D, VectorView<T> E)
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
        // d>0: mu = c - d + d*sqrt(1+|b|^2/d^2)
        //         = c - d*(-|b|^2/d^2) / (1+sqrt(1+|b|^2/d^2)
        //         = c + |b|^2/d / (1+sqrt(1+|b|^2/d^2)
        //
        // d<0: mu = c + |d| - |d|*sqrt(1+|b|^2/d^2)
        //         = c + |d|*(-|b|^2/d^2) / (1+sqrt(1+|b|^2/d^2)
        //         = c - |b|^2/|d| / (1+sqrt(1+|b|^2/d^2)
        //         = c + |b|^2/d / (1+sqrt(1+|b|^2/d^2)
        //
        // mu = c + |b|^2/d / (1 + sqrt(1+|b|^2/d^2))
        const ptrdiff_t N = D.size();
        TMVAssert(E.size() == N-1);
        TMVAssert(N > 1);

        RT a = TMV_NORM(D(N-2)) + (N>2 ? TMV_NORM(E(N-3)) : RT(0));
        RT c = TMV_NORM(D(N-1)) + TMV_NORM(E(N-2));
        T b = D(N-2)*E(N-2);
        RT absb = TMV_ABS(b);
        RT d = (c-a)/2;
        if (d == RT(0)) {
            // This shouldn't happen, but worth checking, since it leads
            // to nan's if we don't check.
            return c+absb;
        } else {
            RT bod = absb/d;
            RT x = absb*bod/(RT(1) + TMV_SQRT(RT(1)+TMV_NORM(bod)));
            return c+x;
        }
    }

    template <class T> 
    static void ReduceBidiagonal22(
        MatrixView<T> U, VectorView<RT> D, VectorView<RT> E, MatrixView<T> Vt)
    {
        // Find Givens rotations which diagonalize 2x2 matrix exactly:
        //
        // [  c1 s1 ] [ d0 e ] [  c2 s2 ] = [  c1 d0   c1 e + s1 d1 ] [  c2 s2 ]
        // [ -s1 c1 ] [ 0 d1 ] [ -s2 c2 ]   [ -s1 d0  -s1 e + c1 d1 ] [ -s2 c2 ]
        // = [ c1 c2 d0 - c1 s2 e - s1 s2 d1    c1 s2 d0 + c1 c2 e + s1 c2 d1 ]
        //   [ -s1 c2 d0 + s1 s2 e - c1 s2 d1  -s1 s2 d0 - s1 c2 e + c1 c2 d1 ]
        //
        // lower left = 0 => s1 (s2 e - c2 d0) = c1 s2 d1
        //                   t1 = t2 d1 / (t2 e - d0)
        //
        // upper right = 0 => s1 c2 d1 = -c1 s2 d0 - c1 c2 e
        //                    t1 = -(t2 d0 + e) / d1
        // 
        // Equate these to get:
        // t2 d1^2 = -(t2 d0 + e)(t2 e - d0) =
        //       -t2^2 d0 e - t2 e^2 + d0^2 t2 + e d0
        // d0 e t2^2 + (d1^2 + e^2 - d0^2) t2 - e d0 = 0
        // t2^2 + ((c-a)/b) t2 - 1 = 0
        // where a = d0^2, b = d0 e, c = d1^2 + e^2 
        // are the elements of B^T B.  (c.f. BidiagonalTrailingEigenvalue)
        //
        // The solution is 
        //    t2 = b / ( d +- sqrt(d^2 + b^2) )
        // where d = (c-a)/2 and we choose the smaller solution.
        //
        // From this, take s2 = b * sign(d), c2 = |d| + sqrt(d^2 + b^2)
        // and then renormalize by sqrt(s2*s2 + c2*c2).
        //
        // Then t1 = -(t2 d0 + e) / d1
        // Do the same to get s1,c1
        // 
        // D(0) <- c1 c2 d0 - c1 s2 e - s1 s2 d1
        //       = c1 c2 ( d0 - t2 e - t1 t2 d1 )
        //       = c1 c2 ( d0 - t2 e - t2 (-t2 d0 - e) )
        //       = c1 c2 ( d0 - t2 e + t2^2 d0 + t2 e )
        //       = c1 c2 d0 (1 + t2^2)
        //       = c1/c2 d0
        // Similarly,
        // D(1) <- c2/c1 d1
        // E(0) <- 0
        //
        // However, these formulae are not sufficiently accurate
        // for the typical case of small b/d.
        //
        // Let b' = b/abs(d)
        //     a' = a/abs(d)
        //     z = sqrt(1+b'^2)
        //
        // t2 = b' / (sign(d) +- z) 
        //
        // if (d > 0): choose + solution.
        //
        //     t2 = b' / (1+z)
        //     t2^2 = b'^2 / (1 + 2z + z^2))
        //     Lemma: x = a/b ==> 1/(1+x) = b/(b+a) = 1 - a/(a+b)
        //     1/(1+t2^2) = 1 - b'^2/(1+b'^2+2z+z^2) = 1 - b'^2/(2z^2+2z)
        //                = 1 - b'^2/(2z(1+z)) = c2^2 = 1-s2^2
        //     s2 = b' / sqrt(2z*(1+z))
        //
        // if (d < 0): choose - solution.
        //
        //     t2 = b' / (-1-z) = -b' / (1+z)
        //     s2 = -b' / sqrt(2z*(1+z))
        //
        // Otherwise, use the same prescription as above.

        TMVAssert(D.size() == 2);
        TMVAssert(E.size() == 1);
        if (U.cptr()) TMVAssert(U.rowsize() == 2);
        if (Vt.cptr()) TMVAssert(Vt.colsize() == 2);

        RT d0 = D(0);
        RT d1 = D(1);
        RT e = E(0);
        dbgcout<<"d0,d1,e = "<<d0<<','<<d1<<','<<e<<std::endl;

        // Rescale to help avoid underflow:
        RT max = TMV_MAX(TMV_ABS(d0),TMV_MAX(TMV_ABS(d1),TMV_ABS(e)));
        d0 /= max;
        d1 /= max;
        e /= max;
        dbgcout<<"d0,d1,e => "<<d0<<','<<d1<<','<<e<<std::endl;

        // If e is small coming in, then this calculation will be 
        // the best we can do to reduce it to zero.  
        // However, if it is largish, then there are some cases where
        // rounding errors keep if from going all the way to zero.
        // So in this case, just do the calculation for what E(0) should
        // be and let ChopSmallElements figure out whether it is small 
        // enough to go the rest of the way to 0.
        //bool exact = (TMV_ABS(e) < 1.e-3);
        // Do it this way, so nan will also trigger an exact answer,
        // which will set e to zero on exit.
        bool exact = !(TMV_ABS(e) > 1.e-3);

        RT d = ((d1-d0)*(d1+d0)+e*e)/RT(2);
        RT absd = TMV_ABS(d);
        dbgcout<<"d,absd = "<<d<<','<<absd<<std::endl;

        RT c1,c2,s1,s2;
        // Usually, d is a largish value compared to the other numbers
        // being calculated, so dividing by absd is a good thing to do.
        // However, if d is small compared to d0^2, d1^2, and e^2, then
        // it usually means that some cancellation happened, and the 
        // normal calculation has significant inaccuracies.  So do the 
        // alternate calculation below which is specialized for small d.
        if (absd > 0.1) {
            dbgcout<<"absd = "<<absd<<std::endl;
            RT b = d0*e/absd;  // This is b' above
            RT z = TMV_SQRT(RT(1)+b*b);
            dbgcout<<"b,z = "<<b<<','<<z<<std::endl;

            s2 = b / TMV_SQRT(RT(2)*z*(z+RT(1))); if (d < 0) s2 = -s2;
            c2 = TMV_SQRT(RT(1)-s2*s2);
            dbgcout<<"s2,c2 = "<<s2<<','<<c2<<std::endl;
            // For small s, improve calculation of c:
            // s^2 = 1-c^2 = (1-c)(1+c)
            // c = 1 - s^2/(1+c)
            if (TMV_ABS(s2) < RT(0.1)) c2 = RT(1) - s2*s2/(RT(1)+c2);
            dbgcout<<"c2 => "<<c2<<std::endl;

            // t1 = t2 d1 / (t2 e - d0) = s2 d1 / (s2 e - d0 c2)
            s1 = s2 * d1;
            c1 = s2 * e - c2 * d0;
            if (TMV_ABS(c1) < 0.1 * TMV_ABS(s1)) {
                dbgcout<<"Small c1 = "<<c1<<std::endl;
                // Then do a calculation that is more accurate for small c1.
                // c1 = s2 e - c2 d0 = e*(s2 - d0/e sqrt(1-s2^2))
                //    = (s2^2 - (d0^2/e^2) (1-s2^2)) /
                //           (s2/e + d0/e^2 sqrt(1-s2^2))
                // The numerator =
                // b^2 / (2z(1+z)) - (d0^2/e^2) + (d0^2/e^2) (b^2 / (2z(1+z)))
                // d0^2 e^2 / (2d^2 z(1+z)) - d0^2/e^2 + d0^4 / (2d^2 z(1+z))
                // d0^2 [ e^2 - 2d^2 z(1+z)/e^2 + d0^2 ] / (2d^2 z(1+z))
                // [] = e^2 + d0^2 - (d1^2-d0^2+e^2) (d z(1+z) / e^2)
                // The key here is that the () term can be ~= 1, so we need to
                // expand it as 1 - (1-()) and cancel the e^2 terms,
                // as |e|~=1 (or often exactly 1 from the normalization).
                // For now, let W = 1-d z(1+z)/e^2.
                // [] = e^2 + d0^2 - (d1^2-d0^2+e^2) (1-W)
                //    = e^2 + d0^2 - d1^2 + d0^2 - e^2 + 2Wd
                //    = 2d0^2 - d1^2 + 2Wd
                // So, numerator =
                // d0^2 (2d0^2 - d1^2 + 2Wd) / (2d^2 z(1+z))
                // And then,
                // c1 = (d0^2 e^2) (2d0^2 - d1^2 + 2Wd)) /
                //         (s2 e + c2 d0) * (2 d^2 z(1+z))
                // 
                // Now, find good expression for W:
                // W = 1 - d z(1+z)/e^2 
                // e^2 W = e^2 - (d1^2-d0^2+e^2) z (1+z)/2
                //       = e^2 - (d1^2-d0^2+e^2) (1-(1-z(1+z)/2))
                //       = e^2 - (d1^2-d0^2+e^2) + 2d (1-z(1+z)/2)
                //       = d0^2-d1^2 + d(2-z-z^2)
                //       = d0^2-d1^2 + d(2-z-(1+b^2))
                //       = d0^2-d1^2 + d(1-z-b^2)
                //       = d0^2-d1^2 - d b^2 + d(1-z^2)/(1+z)
                //       = d0^2-d1^2 - d b^2 + d(-b^2)/(1+z)
                //       = d0^2-d1^2 - d b^2 (1+1/(1+z))
                dbgcout<<"Easy W = "<<1.-d*z*(1.+z)/(e*e)<<std::endl;
                RT W = ((d0-d1)*(d0+d1) - d*b*b*(RT(1)+RT(1)/(RT(1)+z)))/(e*e);
                dbgcout<<"Better W = "<<W<<std::endl;
                RT d0sq = d0*d0/d;
                RT esq = e*e/d;
                c1 = d0sq * esq * (d0*d0+(d0-d1)*(d0+d1)+RT(2)*W*d) /
                    ((s2*e+c2*d0) * (RT(2)*z*(RT(1)+z)));
            }

            RT norm1 = TMV_SQRT(s1*s1 + c1*c1);
            if (c1 < RT(0)) norm1 = -norm1;
            s1 /= norm1;
            c1 /= norm1;
            dbgcout<<"s1,c1 = "<<s1<<','<<c1<<std::endl;
            if (TMV_ABS(s1) < RT(0.1)) c1 = RT(1) - s1*s1/(RT(1)+c1);
            dbgcout<<"c1 => "<<c1<<std::endl;
        } else {
            dbgcout<<"d ~= 0\n";
            // d = (d1^2 + e^2 - d0^2)/2
            // So this means that d0^2 ~= d1^2 + e^2
            // Thus, it is safe to assume that d0 is not small.
            //
            RT b = e/d0; // Now scaling by d0^2 instead of |d|
            d /= d0*d0;
            RT f = d1/d0;
            dbgcout<<"b,d = "<<b<<','<<d<<std::endl;

            // t2 = b sgn(d) / ( |d| + sqrt(d^2 + b^2) )
            s2 = b; if (d < RT(0)) s2 = -s2;
            c2 = absd + TMV_SQRT((d*d + b*b));
            dbgcout<<"s2,c2 = "<<s2<<','<<c2<<std::endl;
            RT norm2 = TMV_SQRT(s2*s2 + c2*c2);
            s2 /= norm2;
            c2 /= norm2;
            dbgcout<<"s2,c2 => "<<s2<<','<<c2<<std::endl;

            // We have two formulae for t1:
            // t1_a = -(t2 + b) / f
            // t1_b = t2 f / (t2 b - 1)
            // For some reason, I'm finding that this branch doesn't
            // always produce suitably accurate answers, but I can't 
            // figure out where the inaccuracy is coming from.  I don't
            // see any place where two largish values nearly cancel.
            // However, we can improve the accuracy by combining the 
            // above two equations to get an equation for t2 in terms of
            // the existing estimate of t2. 
            exact = false;  // And don't assume answer is exact.
            
            // t1 = -(t2 + b) / f
            // t2 = -t1 f - b
            //     = -t2 f^2 / (t2 b - 1) - b
            //     = (-t2 f^2 - t2 b^2 + b) / (t2 b - 1)
            //     = (t2 f^2 + t2 b^2 - b) / (1 - t2 b)
            // t2 (1-t2 b) = t2 (f^2+b^2) - b
            // t2 (f^2+b^2+t2 b - 1) = b
            // This equation is correct for the true t2 = t2~.  
            // But we have a current estimate t2 which is t2~-dt
            // (t2+dt) (f^2+b^2+(t2+dt) b - 1) = b
            // t2 (f^2+b^2+t2b-1) + dt(f^2+b^2+t2b-1) +dt t2b = b
            // dt = [ b - t2(f^2+b^2+t2b-1) ] / [ f^2+b^2+2t2 b-1 ]
            //
            // tan(x+dx) = t + dt
            // tan(x) + dx sec^2(x) = t + dt
            // dx = dt c^2
            // cos(x+dx) = cos(x) - dx sin(x)
            // sin(x+dx) = sin(x) + dx cos(x)
            //
            // So c' = c - s c^2 dt
            //    s' = s + c^3 dt
            RT dt;
            RT t2 = s2/c2;
            dbgcout<<"current t2 = "<<t2<<std::endl;
            dt = ( b - t2*(b*(t2+b)-(1-f)*(1+f)) ) / 
                (b*(2*t2+b)-(1-f)*(1+f));
            dbgcout<<"dt = "<<dt<<std::endl;
            dt *= c2*c2;
            RT ds = c2*dt;
            RT dc = s2*dt;
            dbgcout<<"ds,dc = "<<ds<<','<<dc<<std::endl;
            s2 += ds;
            c2 -= dc;
            dbgcout<<"New s2,c2 = "<<s2<<','<<c2<<std::endl; 

            // Make sure s2^2 + c2^2 is still = 1
            dbgcout<<"s2^2+c2^2 - 1 = "<<s2*s2+c2*c2-RT(1)<<std::endl;
            norm2 = TMV_SQRT(s2*s2 + c2*c2);
            s2 /= norm2;
            c2 /= norm2;
            dbgcout<<"s2,c2 => "<<s2<<','<<c2<<std::endl;

            s1 = s2*f;
            c1 = s2*b-c2;
            dbgcout<<"s1 = "<<-s2<<" + "<<-c2*b<<std::endl;
            dbgcout<<"s1,c1 = "<<s1<<','<<c1<<std::endl;
            RT norm1 = TMV_SQRT(s1*s1 + c1*c1);
            if (c1 < RT(0)) norm1 = -norm1;
            s1 /= norm1;
            c1 /= norm1;
            dbgcout<<"s1,c1  => "<<s1<<','<<c1<<std::endl;
        }
#ifdef XDEBUG
        Matrix<RT> B(2,2); B(0,0) = D(0); B(0,1) = E(0); B(1,0) = RT(0); B(1,1) = D(1);
        Matrix<RT> g1(2,2); g1(0,0) = c1; g1(0,1) = s1; g1(1,0) = -s1; g1(1,1) = c1;
        Matrix<RT> g2(2,2); g2(0,0) = c2; g2(0,1) = s2; g2(1,0) = -s2; g2(1,1) = c2;
        Matrix<RT> S = g1 * B * g2;
        Matrix<T> A(U.cptr()&&Vt.cptr() ? U.colsize() : 0, 
                    U.cptr()&&Vt.cptr() ? Vt.rowsize() : 0);
        if (U.cptr() && Vt.cptr()) A = U * B * Vt;
        dbgcout<<"c1,s1 = "<<c1<<','<<s1<<std::endl;
        dbgcout<<"c2,s2 = "<<c2<<','<<s2<<std::endl;
        dbgcout<<"Initial B = "<<B<<endl;
        dbgcout<<"g1 = "<<g1<<std::endl;
        dbgcout<<"g2 = "<<g2<<std::endl;
        dbgcout<<"S = g1 B g2 = "<<S<<std::endl;
        //if (U.cptr() && Vt.cptr()) dbgcout<<"Initial UBVt = "<<A<<endl;
#endif
        if (exact) {
            dbgcout<<"Use Exact solution for D,E\n";
            D(0) *= c1/c2;
            D(1) *= c2/c1;
            E(0) = RT(0);
        } else {
            dbgcout<<"Use Calculated solution for D,E\n";
            d0 = D(0);
            d1 = D(1);
            e = E(0);
            D(0) = c1*c2*d0 - c1*s2*e - s1*s2*d1;
            D(1) = -s1*s2*d0 - s1*c2*e + c1*c2*d1;
            E(0) = c1*s2*d0 + c1*c2*e + s1*c2*d1;
        }

        if (U.cptr()) {
            Givens<RT> G1(c1,s1);
            G1.mult(U.transpose());
        }
        if (Vt.cptr()) {
            Givens<RT> G2(c2,-s2);
            G2.mult(Vt);
        }
#ifdef XDEBUG
        if (U.cptr() && Vt.cptr()) {
            B.diag() = D;
            B.diag(1) = E;
            Matrix<T> A2 = U * B * Vt;
            dbgcout<<"Done: B = "<<B<<endl;
            dbgcout<<"B-S = "<<B-S<<endl;
            dbgcout<<"Norm(B-S) = "<<Norm(B-S)<<endl;
            //dbgcout<<"UBVt = "<<A2<<endl;
            //dbgcout<<"A2-A = "<<Matrix<T>(A2-A).clip(0.1*MaxAbsElement(A2))<<endl;
            dbgcout<<"Done 22: Norm(A2-A) = "<<Norm(A2-A)<<endl;
            dbgcout<<"Norm(A) = "<<Norm(A)<<std::endl;
            if (!(Norm(A2-A) < THRESH*Norm(A))) {
                cerr<<"ReduceBidiagonal22\n";
                cerr<<"B = "<<B<<endl;
                cerr<<"A = "<<A<<endl;
                cerr<<"S = "<<S<<endl;
                cerr<<"G1 = "<<g1<<endl;
                cerr<<"G2 = "<<g2<<endl;
                cerr<<"A2 = "<<A2<<endl;
                cerr<<"diff = "<<A2-A<<endl;
                cerr<<"diff.maxAbsElement = "<<(A2-A).maxAbsElement()<<std::endl;
                cerr<<"Norm(diff) = "<<Norm(A2-A)<<endl;
                cerr<<"cf. Norm(A) = "<<Norm(A)<<endl;
                cerr<<"THRESH * Norm(A) = "<<THRESH*Norm(A)<<endl;
                abort();
            }
        }
#endif
    }

    template <class T> 
    static void ReduceBidiagonal(
        MatrixView<T> U, VectorView<RT> D, VectorView<RT> E, MatrixView<T> Vt)
    {
#ifdef XDEBUG
        dbgcout<<"Start Reduce Bidiagonal QR:\n";
        dbgcout<<"D = "<<D<<endl;
        dbgcout<<"E = "<<E<<endl;
        Matrix<RT> B(D.size(),D.size(),RT(0));
        Vector<RT> D0 = D;
        Vector<RT> E0 = E;
        B.diag() = D0;
        B.diag(1) = E0;
        const ptrdiff_t M = U.cptr()&&Vt.cptr() ? U.colsize() : 0;
        const ptrdiff_t Nx = U.cptr()&&Vt.cptr() ? Vt.rowsize() : 0;
        Matrix<T> A0(M,Nx);
        if (U.cptr() && Vt.cptr()) {
            A0 = U * B * Vt;
            dbgcout<<"A0 = "<<A0<<endl;
        }
#endif
        // Reduce the superdiagonal elements of Bidiagonal Matrix B 
        // (given by D,E) while maintaining U B Vt. 
        // Note: the input B must be unreduced - ie. all entries are non-zero.
        const ptrdiff_t N = D.size();
        TMVAssert(N>0);
        TMVAssert(E.size() == N-1);
        if (U.cptr()) TMVAssert(U.rowsize() == N);
        if (Vt.cptr()) TMVAssert(Vt.colsize() == N);
        TMVAssert(D.step()==1);
        TMVAssert(E.step()==1);
#ifdef XDEBUG
        TMVAssert(D.minAbs2Element() > RT(0));
        TMVAssert(E.minAbs2Element() > RT(0));
#endif

        if (N == 1) return;
        if (N == 2) return ReduceBidiagonal22(U,D,E,Vt);

        // The reduction is based on the QR algorithm to diagonalize the
        // unreduced symmetric tridiagonal matrix T = BtB
        // The basic idea is as follows:
        // (see Golub and van Loan, chapter 8 for a derivation)
        //
        // if T is a symmetric tridiagonal matrix
        // and mu is (approximately) an eigenvalue of T
        // and the QR decomposition of T - mu I = Vt R
        // Then, T' = R Vt + mu I will be tridiagonal with the last 
        // subdiagonal element small.
        // (Note: T' = V (T-muI) Vt + muI = V T Vt.)
        //
        // Wilkinson (1968) suggested that a good choice for mu is
        // the eigenvalue of the trailing 2x2 block of T that is 
        // closer to the trailing diagonal element.
        //
        // Rather than explicitly forming T = BtB and doing this
        // procedure, Golub and van Load show that it can be done
        // in place.
        // If T' = V T Vt, 
        // then Bt'B' = V Bt B Vt
        // B' = U B Vt for some U
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
        // For each Givens rotation we use, we also multiply U or Vt by the
        // adjoint to maintain the constancy of U B Vt.
        //
        // At the end of this procedure, E(N-1) should be smaller than it was.
        // Note: This procedure works exactly if N=2.
        RT* Di = D.ptr();
        RT* Ei = E.ptr();

        RT mu = BidiagonalTrailingEigenValue(D,E);
        dbgcout<<"mu = "<<mu<<std::endl;
        RT y = TMV_NORM(*Di) - mu;  // = T00 - mu
        dbgcout<<"y = "<<y<<std::endl;
        RT x = TMV_CONJ(*Di)*(*Ei);  // = T10
        dbgcout<<"x = "<<x<<std::endl;
        Givens<RT> G = GivensRotate(y,x);
        dbgcout<<"Rotatedi y,x => "<<y<<','<<x<<std::endl;
        for(ptrdiff_t i=1;i<N;++i) {
            G.mult(*Di,*Ei);
            dbgcout<<"D,E -> "<<*Di<<','<<*Ei<<std::endl;
            if (Vt.cptr()) G.conjMult(Vt.rowPair(i-1,i));
            TMVAssert(x==RT(0));
            G.mult(x,*(++Di)); // x = B(i,i-1)
            dbgcout<<"x,D -> "<<x<<','<<*Di<<std::endl;
            G = GivensRotate(*(Di-1),x);
            dbgcout<<"Rotated D,x => "<<*(Di-1)<<','<<x<<std::endl;
            G.mult(*Ei,*Di);
            dbgcout<<"E,D -> "<<*Ei<<','<<*Di<<std::endl;
            if (U.cptr()) G.conjMult(U.colPair(i-1,i).transpose());
            if (i < N-1) {
                TMVAssert(x==RT(0));
                G.mult(x,*(++Ei)); // x = B(i-1,i+1)
                dbgcout<<"x,E -> "<<i<<','<<*Ei<<std::endl;
                G = GivensRotate(*(Ei-1),x);
                dbgcout<<"Rotated E,x => "<<*(Ei-1)<<','<<x<<std::endl;
            }
        }
#ifdef XDEBUG
        if (U.cptr() && Vt.cptr()) {
            B.diag() = D;
            B.diag(1) = E;
            Matrix<T> AA = U * B * Vt;
            if (!(Norm(A0-AA) < THRESH*Norm(A0))) {
                cerr<<"ReduceBidiagonal: \n";
                cerr<<"input D = "<<D0<<endl;
                cerr<<"input E = "<<E0<<endl;
                cerr<<"UBVt = "<<A0<<endl;
                cerr<<"UBVt => "<<AA<<endl;
                cerr<<"U = "<<U<<endl;
                cerr<<"D = "<<D<<endl;
                cerr<<"E = "<<E<<endl;
                cerr<<"Vt = "<<Vt<<endl;
                abort();
            }
        }
        if (Norm(D-D0) + Norm(E-E0) == RT(0)) {
            cerr<<"ReduceBidiagonal didn't reduce:\n";
            cerr<<"input D = "<<D0<<endl;
            cerr<<"input E = "<<E0<<endl;
            cerr<<"output D = "<<D<<endl;
            cerr<<"output E = "<<E<<endl;
            cerr<<"mu = "<<mu<<endl;
            abort();
        }
#endif
    }

    template <class T> 
    void SV_DecomposeFromBidiagonal_QR(
        MatrixView<T> U, VectorView<RT> D, VectorView<RT> E, MatrixView<T> Vt)
    {
#ifdef XDEBUG
        dbgcout<<"Start Decompose from Bidiagonal QR:\n";
        //if (U.cptr()) dbgcout<<"U = "<<U<<endl;
        //if (Vt.cptr()) dbgcout<<"Vt = "<<Vt<<endl;
        dbgcout<<"D = "<<D<<endl;
        dbgcout<<"E = "<<E<<endl;
        Matrix<RT> B(D.size(),D.size(),RT(0));
        B.diag() = D;
        B.diag(1) = E;
        const ptrdiff_t M = U.cptr()&&Vt.cptr() ? U.colsize() : 0;
        const ptrdiff_t Nx = U.cptr()&&Vt.cptr() ? Vt.rowsize() : 0;
        Matrix<T> A0(M,Nx);
        if (U.cptr() && Vt.cptr()) {
            A0 = U * B * Vt;
            dbgcout<<"A0 = "<<A0<<endl;
        }
#endif

        const ptrdiff_t N = D.size();
        if (N <= 1) return;

        TMVAssert(D.size() == E.size()+1);
        if (U.cptr()) TMVAssert(U.rowsize() == D.size());
        if (Vt.cptr()) TMVAssert(Vt.colsize() == D.size());
#ifdef XDEBUG
        TMVAssert(D.minAbs2Element() > RT(0));
        TMVAssert(E.minAbs2Element() > RT(0));
#endif

        // We successively reduce the superdiagonal of B (E) to 0
        // using a sequence of Givens rotations. 
        // The reduction procedure tends to push the values up and left, so it 
        // makes sense to start at the lower right and work back up the matrix.
        // We also set to zero any very small values based on machine precision.
        // Loop invariant: all E(i) with i>=q are 0.
        // Initially q = N-1. (ie. All E(i) are potentially non-zero.)
        // When q = 0, we are done.

        for(ptrdiff_t q=N-1; q>0; ) {
            if (E(q-1) == T(0)) --q;
            else {
                ptrdiff_t p=q-1;
                while (p > 0 && (E(p-1) != T(0))) --p;
                // Set p such that E(p-1) = 0 and all E(i) with p<=i<q are 
                // non-zero.

                if (U.cptr()) {
                    if (Vt.cptr()) ReduceBidiagonal<T>(
                        U.colRange(p,q+1),D.subVector(p,q+1),
                        E.subVector(p,q),Vt.rowRange(p,q+1));
                    else ReduceBidiagonal<T>(
                        U.colRange(p,q+1),D.subVector(p,q+1),
                        E.subVector(p,q),Vt);
                } else {
                    if (Vt.cptr()) ReduceBidiagonal<T>(
                        U,D.subVector(p,q+1),
                        E.subVector(p,q),Vt.rowRange(p,q+1));
                    else ReduceBidiagonal<T>(
                        U,D.subVector(p,q+1),E.subVector(p,q),Vt);
                }

                bool newzeroD = false;
#ifdef XDEBUG
                dbgcout<<"Before Chop: \n";
                dbgcout<<"D = "<<D.subVector(p,q+1)<<std::endl;
                dbgcout<<"E = "<<E.subVector(p,q)<<std::endl;
#endif
                BidiagonalChopSmallElements(
                    D.subVector(p,q+1),E.subVector(p,q),&newzeroD);
#ifdef XDEBUG
                dbgcout<<"After Chop: \n";
                dbgcout<<"D = "<<D.subVector(p,q+1)<<std::endl;
                dbgcout<<"E = "<<E.subVector(p,q)<<std::endl;
                dbgcout<<"newzero in D? "<<newzeroD<<std::endl;
#endif
                // Check that we haven't introduced new 0's in the D vector.
                // If we have, we need to go back to the original calling 
                // routine, which deals with them.
                if (newzeroD) {
                    if (U.cptr()) {
                        if (Vt.cptr()) DoSVDecomposeFromBidiagonal<T>(
                            U.colRange(p,q+1),D.subVector(p,q+1),
                            E.subVector(p,q),Vt.rowRange(p,q+1),false,false);
                        else DoSVDecomposeFromBidiagonal<T>(
                            U.colRange(p,q+1),D.subVector(p,q+1),
                            E.subVector(p,q),Vt,false,false);
                    } else {
                        if (Vt.cptr()) DoSVDecomposeFromBidiagonal<T>(
                            U,D.subVector(p,q+1),
                            E.subVector(p,q),Vt.rowRange(p,q+1),false,false);
                        else DoSVDecomposeFromBidiagonal<T>(
                            U,D.subVector(p,q+1),
                            E.subVector(p,q),Vt,false,false);
                    }
                    q = p;
                }
            }
        }
#ifdef XDEBUG
        if (U.cptr() && Vt.cptr()) {
            Matrix<T> AA = U * DiagMatrixViewOf(D) * Vt;
            dbgcout<<"Done QR Norm(A0-AA) = "<<Norm(A0-AA)<<endl;
            dbgcout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<endl;
            dbgcout<<"Norm(VtV-1) = "<<Norm(Vt*Vt.adjoint()-T(1))<<endl;
            dbgcout<<"Norm(VVt-1) = "<<Norm(Vt.adjoint()*Vt-T(1))<<endl;
            dbgcout<<"U = "<<U<<std::endl;
            dbgcout<<"S = "<<D<<std::endl;
            dbgcout<<"Vt = "<<Vt<<std::endl;
            if (!(Norm(A0-AA) < THRESH*Norm(A0))) {
                cerr<<"SV_DecomposeFromBidiagonal QR: \n";
                cerr<<"UBVt = "<<A0<<endl;
                cerr<<"USVt = "<<AA<<endl;
                //cerr<<"U = "<<U<<endl;
                //cerr<<"Vt = "<<Vt<<endl;
                cerr<<"input B = "<<B<<endl;
                cerr<<"S = "<<D<<endl;
                dbgcout<<"Norm(UBVt-USVt) = "<<Norm(A0-AA)<<endl;
                dbgcout<<"cf "<<THRESH<<"*Norm(UBVt) = "<<THRESH*Norm(A0)<<endl;
                abort();
            }
        }
#endif
    }

#undef RT

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_SVDecompose_QR.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


