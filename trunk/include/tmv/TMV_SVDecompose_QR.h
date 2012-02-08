

#ifndef TMV_SVDecompose_QR_H
#define TMV_SVDecompose_QR_H

#include "TMV_SVDecompose.h"

#ifdef XDEBUG_SVD
#define THRESH 1.e-4
#include "TMV_ProdMM.h"
#include "TMV_MultMM.h"
#include "TMV_SumMM.h"
#include "TMV_SumVV.h"
#include "TMV_SumMX.h"
#include "TMV_AddMM.h"
#include "TMV_Norm.h"
#include "TMV_NormM.h"
#endif

#ifdef PRINTALGO_SVD
#include <iostream>
#include "TMV_MatrixIO.h"
#include "TMV_VectorIO.h"
#define TMV_USE_WRITER
#endif

#define dbgcout SafeWriter()
#include "TMV_SafeWriter.h"

namespace tmv {

    // Defined in TMV_SVDecompose_QR.cpp
    template <class Tu, class RT, class Tv>
    void InstSV_DecomposeFromBidiagonal_QR(
        MatrixView<Tu> U, VectorView<RT> D, VectorView<RT> E,
        MatrixView<Tv> V, bool UisI, bool VisI);

    // Return the Wilkinson choice for an eigenvalue of T = BtB, namely
    // the eigenvalue of the trailing 2x2 block of T which is closer
    // to the last diagonal element of T.
    template <class Vd, class Ve>
    static typename Vd::real_type BidiagonalTrailingEigenValue(
        BaseVector_Mutable<Vd>& D, BaseVector_Mutable<Ve>& E)
    {
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
        typedef typename Vd::value_type T;
        typedef typename Vd::real_type RT;

        RT a = TMV_NORM(D.cref(N-2)) + (N>2 ? TMV_NORM(E.cref(N-3)) : RT(0));
        RT c = TMV_NORM(D.cref(N-1)) + TMV_NORM(E.cref(N-2));
        T b = D.cref(N-2)*E.cref(N-2);
        RT absb = TMV_ABS(b);
        RT d = (c-a)/2;
        if (d == RT(0)) {
            // This shouldn't happen, but worth checking, since it leads
            // to nan's if it does happen.
            return c+absb;
        } else {
            RT bod = absb/d;
            RT x = absb*bod/(RT(1) + TMV_SQRT(RT(1)+TMV_NORM(bod)));
            return c+x;
        }
    }
                        
    template <int algo, ptrdiff_t cs, ptrdiff_t rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_QR_Helper;

    // algo 0: Trivial, nothing to do (M == 0, or N == 0)
    template <ptrdiff_t cs, ptrdiff_t rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_QR_Helper<0,cs,rs,Mu,Vd,Ve,Mv>
    { static TMV_INLINE void call(Mu& , Vd& , Ve& , Mv& , bool , bool ) {} };

    // algo 11: Normal algorithm
    template <ptrdiff_t cs, ptrdiff_t rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_QR_Helper<11,cs,rs,Mu,Vd,Ve,Mv>
    {
        static void call(Mu& U, Vd& D, Ve& E, Mv& V, bool UisI, bool VisI)
        {
            typedef typename Vd::value_type RT;
            TMVStaticAssert(Traits<RT>::isreal);

            const ptrdiff_t N = rs==Unknown ? D.size() : rs;
            if (N <= 1) return;
#ifdef PRINTALGO_SVD
            const ptrdiff_t M = cs==Unknown ? U.colsize() : cs;
            std::cout<<"SVDecomposeFromBidiagonal algo 11: M,N,cs,rs = "<<
                M<<','<<N<<','<<cs<<','<<rs<<std::endl;
#endif
            const ptrdiff_t xx = Unknown;
            typedef typename Mu::colrange_type Mus;
            typedef typename Vd::subvector_type Vds;
            typedef typename Ve::subvector_type Ves;
            typedef typename Mv::rowrange_type Mvs;

            // We successively reduce the superdiagonal of B (E) to 0
            // using a sequence of Givens rotations. 
            // The reduction procedure tends to push the values up and left, 
            // so it makes sense to start at the lower right and work back 
            // up the matrix.
            // We also set to zero any very small values based on machine 
            // precision.
            // Loop invariant: all E(i) with i>=q are 0.
            // Initially q = N-1. (ie. All E(i) are potentially non-zero.)
            // When q = 0, we are done.

            for(ptrdiff_t q=N-1; q>0; ) {
                if (E.cref(q-1) == RT(0)) --q;
                else {
                    ptrdiff_t p=q-1;
                    while (p > 0 && (E(p-1) != RT(0))) --p;
                    // Set p such that E(p-1) = 0 and all E(i) with p<=i<q are 
                    // non-zero.

                    Mus U1 = U.cptr() ? U.colRange(p,q+1) : U.colRange(0,0);
                    Vds D1 = D.subVector(p,q+1);
                    Ves E1 = E.subVector(p,q);
                    Mvs V1 = V.cptr() ? V.rowRange(p,q+1) : V.rowRange(0,0);
                    SVDecomposeFromBidiagonal_QR_Helper<
                        21,xx,xx,Mus,Vds,Ves,Mvs>::reduce(U1,D1,E1,V1);

                    bool newzeroD = false;
#ifdef XDEBUG_SVD
                    dbgcout<<"Before Chop: \n";
                    dbgcout<<"D = "<<D1<<std::endl;
                    dbgcout<<"E = "<<E1<<std::endl;
#endif
                    BidiagonalChopSmallElements(D1,E1,newzeroD);
#ifdef XDEBUG_SVD
                    dbgcout<<"After Chop: \n";
                    dbgcout<<"D = "<<D1<<std::endl;
                    dbgcout<<"E = "<<E1<<std::endl;
                    dbgcout<<"newzero in D? "<<newzeroD<<std::endl;
#endif
                    // Check that we haven't introduced new 0's in the D vector.
                    // If we have, we need to go back to the original calling 
                    // routine, which deals with them.
                    if (newzeroD) {
                        SVDecomposeFromBidiagonal_Helper<
                            31,xx,xx,Mus,Vds,Ves,Mvs>::call(
                                U1,D1,E1,V1,false,false);
                        q = p;
                    }
                }
            }
        }
    };

    // algo 21: Try to reduce the last elemet of E to zero.
    template <ptrdiff_t cs, ptrdiff_t rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_QR_Helper<21,cs,rs,Mu,Vd,Ve,Mv>
    {
        static void reduce(Mu& U, Vd& D, Ve& E, Mv& V)
        {
            typedef typename Mu::real_type RT;
            const ptrdiff_t N = rs==Unknown ? D.size() : rs;
#ifdef PRINTALGO_SVD
            const ptrdiff_t M = cs==Unknown ? U.colsize() : cs;
            std::cout<<"SVDecomposeFromBidiagonal algo 21: M,N,cs,rs = "<<
                M<<','<<N<<','<<cs<<','<<rs<<std::endl;
#endif
            if (N <= 1) return;
            if (N == 2) 
                return SVDecomposeFromBidiagonal_QR_Helper<
                    22,cs,rs,Mu,Vd,Ve,Mv>::reduce(U,D,E,V);
#ifdef XDEBUG_SVD
            TMVAssert(D.minAbsElement() > RT(0));
            TMVAssert(E.minAbsElement() > RT(0));
            typedef typename Traits2<typename Mu::value_type, typename Mv::value_type>::type T;
            dbgcout<<"Start Reduce Bidiagonal QR:\n";
            dbgcout<<"U = "<<U<<std::endl;
            dbgcout<<"D = "<<D<<std::endl;
            dbgcout<<"E = "<<E<<std::endl;
            dbgcout<<"V = "<<V<<std::endl;
            Matrix<RT> B(N,N,RT(0));
            Vector<RT> D0 = D;
            Vector<RT> E0 = E;
            B.diag() = D;
            B.diag(1) = E;
            const ptrdiff_t M1 = U.cptr() && V.cptr() ? U.colsize() : 0;
            const ptrdiff_t N1 = U.cptr() && V.cptr() ? V.rowsize() : 0;
            Matrix<T> A0(M1,N1);
            if (U.cptr() && V.cptr()) A0 = U * B * V;
            //dbgcout<<"A0 = "<<A0<<std::endl;
#endif
            // Reduce the superdiagonal elements of Bidiagonal Matrix B 
            // (given by D,E) while maintaining U B V. 
            // Note: the input B must be unreduced - ie. all entries are 
            // non-zero.

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
            // For each Givens rotation we use, we also multiply U or V by the
            // adjoint to maintain the constancy of U B V.
            //
            // At the end of this procedure, E(N-1) should be smaller than it was.
            // Note: This procedure should work exactly if N=2.
            // (But it sometimes doesn't.  Hence the special 2x2 
            // algorithm below...)
            
            typedef typename Mu::colpair_type::transpose_type Mucpt;
            typedef typename Mv::rowpair_type Mvrp;
            typedef typename Vd::iterator IT1;
            typedef typename Ve::iterator IT2;
            IT1 Di = D.begin();
            IT2 Ei = E.begin();
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
                if (V.cptr()) {
                    Mvrp Vrp = V.rowPair(i-1,i);
                    G.conjMult(Vrp);
                }
                TMVAssert(x==RT(0));
                G.mult(x,*(++Di)); // x = B(i,i-1)
                dbgcout<<"x,D -> "<<x<<','<<*Di<<std::endl;
                G = GivensRotate(*(Di-1),x);
                dbgcout<<"Rotatedi D,x => "<<*(Di-1)<<','<<x<<std::endl;
                G.mult(*Ei,*Di);
                dbgcout<<"E,D -> "<<*Ei<<','<<*Di<<std::endl;
                if (U.cptr()) {
                    Mucpt Ucpt = U.colPair(i-1,i).transpose();
                    G.conjMult(Ucpt);
                }
                if (i < N-1) {
                    TMVAssert(x==RT(0));
                    G.mult(x,*(++Ei)); // x = B(i-1,i+1)
                    dbgcout<<"x,E -> "<<i<<','<<*Ei<<std::endl;
                    G = GivensRotate(*(Ei-1),x);
                    dbgcout<<"Rotatedi E,x => "<<*(Ei-1)<<','<<x<<std::endl;
                }
            }
#ifdef XDEBUG_SVD
            if (U.cptr() && V.cptr()) {
                B.diag() = D;
                B.diag(1) = E;
                Matrix<T> AA = U * B * V;
                if (!(Norm(A0-AA) <= THRESH*Norm(A0))) {
                    std::cerr<<"ReduceBidiagonal: \n";
                    std::cerr<<"input D = "<<D0<<std::endl;
                    std::cerr<<"input E = "<<E0<<std::endl;
                    std::cerr<<"UBV = "<<A0<<std::endl;
                    std::cerr<<"UBV => "<<AA<<std::endl;
                    std::cerr<<"diff = "<<A0-AA<<std::endl;
                    std::cerr<<"U = "<<U<<std::endl;
                    std::cerr<<"D = "<<D<<std::endl;
                    std::cerr<<"E = "<<E<<std::endl;
                    std::cerr<<"V = "<<V<<std::endl;
                    std::cerr<<"Norm(diff) = "<<Norm(A0-AA)<<std::endl;
                    std::cerr<<"cf. Norm(A) = "<<Norm(A0)<<std::endl;
                    std::cerr<<"THRESH * Norm(A) = "<<THRESH*Norm(A0)<<std::endl;
                    abort();
                }
            }
            if (Norm(D-D0) + Norm(E-E0) == RT(0)) {
                std::cerr<<"ReduceBidiagonal didn't reduce:\n";
                std::cerr<<"input D = "<<D0<<std::endl;
                std::cerr<<"input E = "<<E0<<std::endl;
                std::cerr<<"output D = "<<D<<std::endl;
                std::cerr<<"output E = "<<E<<std::endl;
                std::cerr<<"mu = "<<mu<<std::endl;
                abort();
            }
#endif
        }
    };

    // algo 22: Special case for N = 2
    // Occasionally, with the regular algorithm, rounding errors can
    // lead to an infinite loop.  So this version tries to be more accurate
    // for the 2x2 case, which is mathematically exact.  But we do
    // some special things to make sure rounding errors don't kill
    // the reduction.
    template <ptrdiff_t cs, ptrdiff_t rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_QR_Helper<22,cs,rs,Mu,Vd,Ve,Mv>
    {
        static void reduce(Mu& U, Vd& D, Ve& E, Mv& V)
        {
            typedef typename Mu::real_type RT;
#ifdef PRINTALGO_SVD
            const ptrdiff_t M = cs==Unknown ? U.colsize() : cs;
            const ptrdiff_t N = rs==Unknown ? D.size() : rs;
            std::cout<<"SVDecomposeFromBidiagonal algo 22: M,N,cs,rs = "<<
                M<<','<<N<<','<<cs<<','<<rs<<std::endl;
#endif
            // Find Givens rotations which diagonalize 2x2 matrix exactly:
            //
            // [ c1 s1] [d0 e] [ c2 s2] = [ c1 d0   c1 e + s1 d1] [ c2 s2]
            // [-s1 c1] [0 d1] [-s2 c2]   [-s1 d0  -s1 e + c1 d1] [-s2 c2]
            // = [c1c2d0 - c1s2e - s1s2d1    c1s2d0 + c1c2e + s1c2d1 ]
            //   [-s1c2d0 + s1s2e - c1s2d1  -s1s2d0 - s1c2e + c1c2d1 ]
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
            if (V.cptr()) TMVAssert(V.colsize() == 2);

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
            // Do it this way (with !(>) rather than <), so nan will also 
            // trigger an exact answer, which will set e to zero on exit.
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
                    // b^2/(2z(1+z)) - (d0^2/e^2) + (d0^2/e^2) (b^2/(2z(1+z)))
                    // d0^2 e^2/(2d^2 z(1+z)) - d0^2/e^2 + d0^4/(2d^2 z(1+z))
                    // d0^2 [ e^2 - 2d^2 z(1+z)/e^2 + d0^2 ] / (2d^2 z(1+z))
                    // [] = e^2 + d0^2 - (d1^2-d0^2+e^2) (d z(1+z)/e^2)
                    // The key here is that the () term can be ~= 1, so we 
                    // need to expand it as 1-(1-()) and cancel the e^2 terms,
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
                    RT W = ((d0-d1)*(d0+d1) -
                            d*b*b*(RT(1)+RT(1)/(RT(1)+z)))/(e*e);
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
#ifdef XDEBUG_SVD
            typedef typename Traits2<typename Mu::value_type, typename Mv::value_type>::type T;
            Matrix<RT> B(2,2); 
            B(0,0) = D(0); B(0,1) = E(0); B(1,0) = RT(0); B(1,1) = D(1);
            Matrix<RT> g1(2,2); 
            g1(0,0) = c1; g1(0,1) = s1; g1(1,0) = -s1; g1(1,1) = c1;
            Matrix<RT> g2(2,2); 
            g2(0,0) = c2; g2(0,1) = s2; g2(1,0) = -s2; g2(1,1) = c2;
            Matrix<RT> S = g1 * B * g2;
            const ptrdiff_t M1 = U.cptr() && V.cptr() ? U.colsize() : 0;
            const ptrdiff_t N1 = U.cptr() && V.cptr() ? V.rowsize() : 0;
            Matrix<T> A(M1,N1);
            if (U.cptr() && V.cptr()) A = U * B * V;
            dbgcout<<"c1,s1 = "<<c1<<','<<s1<<std::endl;
            dbgcout<<"c2,s2 = "<<c2<<','<<s2<<std::endl;
            dbgcout<<"Initial B = "<<B<<std::endl;
            dbgcout<<"g1 = "<<g1<<std::endl;
            dbgcout<<"g2 = "<<g2<<std::endl;
            dbgcout<<"S = g1 B g2 = "<<S<<std::endl;
            //if (U && V) dbgcout<<"Initial UBV = "<<A<<std::endl;
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
                typename Mu::transpose_type Ut = U.transpose();
                G1.mult(Ut);
            }
            if (V.cptr()) {
                Givens<RT> G2(c2,-s2);
                G2.mult(V);
            }
#ifdef XDEBUG_SVD
            if (U.cptr() && V.cptr()) {
                B.diag() = D;
                B.diag(1) = E;
                Matrix<T> A2 = U * B * V;
                dbgcout<<"Done: B = "<<B<<std::endl;
                dbgcout<<"B-S = "<<B-S<<std::endl;
                dbgcout<<"Norm(B-S) = "<<Norm(B-S)<<std::endl;
                //dbgcout<<"UBV = "<<A2<<std::endl;
                //dbgcout<<"A2-A = "<<Matrix<T>(A2-A).clip(0.1*MaxAbsElement(A2))<<std::endl;
                dbgcout<<"Done 22: Norm(A2-A) = "<<Norm(A2-A)<<std::endl;
                dbgcout<<"Norm(A) = "<<Norm(A)<<std::endl;
                if (!(Norm(A2-A) <= THRESH*Norm(A))) {
                    std::cerr<<"ReduceBidiagonal22\n";
                    std::cerr<<"B = "<<B<<std::endl;
                    std::cerr<<"A = "<<A<<std::endl;
                    std::cerr<<"S = "<<S<<std::endl;
                    std::cerr<<"G1 = "<<g1<<std::endl;
                    std::cerr<<"G2 = "<<g2<<std::endl;
                    std::cerr<<"A2 = "<<A2<<std::endl;
                    std::cerr<<"diff = "<<A2-A<<std::endl;
                    std::cerr<<"diff.maxAbsElement = "<<(A2-A).maxAbsElement()<<std::endl;
                    std::cerr<<"Norm(diff) = "<<Norm(A2-A)<<std::endl;
                    std::cerr<<"cf. Norm(A) = "<<Norm(A)<<std::endl;
                    std::cerr<<"THRESH * Norm(A) = "<<THRESH*Norm(A)<<std::endl;
                    abort();
                }
            }
#endif

        }
    };

    // algo 90: call InstSV_DecomposeFromBidiagonal_QR
    template <ptrdiff_t cs, ptrdiff_t rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_QR_Helper<90,cs,rs,Mu,Vd,Ve,Mv>
    {
        static TMV_INLINE void call(
            Mu& U, Vd& D, Ve& E, Mv& V, bool UisI, bool VisI)
        { 
            InstSV_DecomposeFromBidiagonal_QR(
                U.xView(),D.xView(),E.xView(),V.xView(),UisI,VisI); 
        }
    };

    // algo -3: Determine which algorithm to use
    template <ptrdiff_t cs, ptrdiff_t rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_QR_Helper<-3,cs,rs,Mu,Vd,Ve,Mv>
    {
        static TMV_INLINE void call(
            Mu& U, Vd& D, Ve& E, Mv& V, bool UisI, bool VisI)
        {
            const int algo = 11;
#ifdef PRINTALGO_SVD
            std::cout<<"Inline SVDecomposeFromBidiagonal: \n";
            std::cout<<"U = "<<TMV_Text(U)<<std::endl;
            std::cout<<"D = "<<TMV_Text(D)<<std::endl;
            std::cout<<"E = "<<TMV_Text(E)<<std::endl;
            std::cout<<"V = "<<TMV_Text(V)<<std::endl;
            std::cout<<"cs = "<<cs<<"  rs = "<<rs<<std::endl;
            std::cout<<"sizes = "<<U.colsize()<<"  "<<U.rowsize()<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
#ifdef XDEBUG_SVD
            typedef typename Traits2<typename Mu::value_type, typename Mv::value_type>::type T;
            typedef typename Mu::real_type RT;

            const ptrdiff_t N = rs==Unknown ? D.size() : rs;
            dbgcout<<"Start Decompose from Bidiag:\n";
            if (U.cptr()) dbgcout<<"U = "<<TMV_Text(U)<<std::endl;
            if (V.cptr()) dbgcout<<"V = "<<TMV_Text(V)<<std::endl;
            dbgcout<<"D = "<<TMV_Text(D)<<"  step "<<D.step()<<"  "<<D<<std::endl;
            dbgcout<<"E = "<<TMV_Text(E)<<"  step "<<E.step()<<"  "<<E<<std::endl;
            //if (U.cptr()) dbgcout<<"U = "<<U<<std::endl;
            //if (V.cptr()) dbgcout<<"V = "<<V<<std::endl;

            dbgcout<<"UisI, VisI = "<<UisI<<"  "<<VisI<<std::endl;
            Matrix<RT> B(N,N,RT(0));
            B.diag() = D;
            B.diag(1) = E;
            const ptrdiff_t M1 = U.cptr() && V.cptr() ? U.colsize() : 0;
            const ptrdiff_t N1 = U.cptr() && V.cptr() ? V.rowsize() : 0;
            Matrix<T> A0(M1,N1);
            if (U.cptr() && V.cptr()) A0 = U * B * V;
            //dbgcout<<"A0 = "<<A0<<std::endl;
#endif
            SVDecomposeFromBidiagonal_QR_Helper<algo,cs,rs,Mu,Vd,Ve,Mv>::call( 
                U,D,E,V,UisI,VisI);
            //std::cout<<"U = "<<U<<std::endl;
            //std::cout<<"S = "<<S<<std::endl;
            //std::cout<<"V = "<<V<<std::endl;
#ifdef XDEBUG_SVD
            if (U.cptr() && V.cptr()) {
                Matrix<T> AA = U * DiagMatrixViewOf(D) * V;
                dbgcout<<"Done QR Norm(A0-AA) = "<<Norm(A0-AA)<<std::endl;
                dbgcout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<std::endl;
                dbgcout<<"Norm(VtV-1) = "<<Norm(V.adjoint()*V-T(1))<<std::endl;
                dbgcout<<"Norm(VVt-1) = "<<Norm(V*V.adjoint()-T(1))<<std::endl;
                dbgcout<<"U = "<<U<<std::endl;
                dbgcout<<"S = "<<D<<std::endl;
                dbgcout<<"V = "<<V<<std::endl;
                if (!(Norm(A0-AA) <= THRESH*Norm(A0))) {
                    std::cerr<<"SV_DecomposeFromBidiagonal QR: \n";
                    std::cerr<<"UBV = "<<A0<<std::endl;
                    std::cerr<<"USV = "<<AA<<std::endl;
                    //std::cerr<<"U = "<<U<<std::endl;
                    //std::cerr<<"V = "<<V<<std::endl;
                    std::cerr<<"input B = "<<B<<std::endl;
                    std::cerr<<"S = "<<D<<std::endl;
                    dbgcout<<"Norm(UBV-USV) = "<<Norm(A0-AA)<<std::endl;
                    dbgcout<<"cf "<<THRESH<<"*Norm(UBV) = "<<THRESH*Norm(A0)<<std::endl;
                    abort();
                }
            }
#endif
        }
    };

    // algo -2: Check for inst
    template <ptrdiff_t cs, ptrdiff_t rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_QR_Helper<-2,cs,rs,Mu,Vd,Ve,Mv>
    {
        static TMV_INLINE void call(
            Mu& U, Vd& D, Ve& E, Mv& V, bool UisI, bool VisI)
        {
            typedef typename Mu::value_type T;
            const bool inst = 
                (cs == Unknown || cs > 16) &&
                (rs == Unknown || rs > 16) &&
                Traits<T>::isinst;
            const int algo = 
                inst ? 90 :
                -3;
            SVDecomposeFromBidiagonal_QR_Helper<algo,cs,rs,Mu,Vd,Ve,Mv>::call(
                U,D,E,V,UisI,VisI);
        }
    };

    template <ptrdiff_t cs, ptrdiff_t rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_QR_Helper<-1,cs,rs,Mu,Vd,Ve,Mv>
    {
        static TMV_INLINE void call(
            Mu& U, Vd& D, Ve& E, Mv& V, bool UisI, bool VisI)
        {
            SVDecomposeFromBidiagonal_QR_Helper<-2,cs,rs,Mu,Vd,Ve,Mv>::call(
                U,D,E,V,UisI,VisI); 
        }
    };

    template <class Mu, class Vd, class Ve, class Mv>
    inline void InlineSV_DecomposeFromBidiagonal_QR(
        BaseMatrix_Rec_Mutable<Mu>& U,
        BaseVector_Mutable<Vd>& D, BaseVector_Mutable<Ve>& E, 
        BaseMatrix_Rec_Mutable<Mv>& V, bool UisI, bool VisI)
    {
        TMVStaticAssert((Traits2<
                         typename Mu::value_type,
                         typename Vd::value_type>::samebase));
        TMVStaticAssert((Traits2<
                         typename Vd::value_type,
                         typename Ve::value_type>::sametype));
        TMVStaticAssert((Traits2<
                         typename Mu::value_type,
                         typename Mv::value_type>::samebase));
        TMVStaticAssert(Traits<typename Vd::value_type>::isreal);
        TMVStaticAssert(!Mu::_conj);
        TMVStaticAssert(!Mv::_conj);
        TMVStaticAssert((Sizes<Mu::_rowsize,Vd::_size>::same));
        TMVStaticAssert((Sizes<Vd::_size,IntTraits<Ve::_size>::Sp1>::same));
        TMVStaticAssert((Sizes<Mu::_rowsize,Mv::_colsize>::same));
        TMVStaticAssert((Sizes<Mu::_rowsize,Mv::_rowsize>::same));
        TMVAssert(D.size() == E.size()+1);
        if (U.cptr()) {
            TMVAssert(U.colsize() >= U.rowsize());
            TMVAssert(U.rowsize() == D.size());
        }
        if (V.cptr()) {
            TMVAssert(V.colsize() == D.size());
            TMVAssert(V.rowsize() >= V.colsize());
        }
        const ptrdiff_t cs = Mu::_colsize;
        const ptrdiff_t rs1 = Sizes<Mu::_rowsize,Vd::_size>::size;
        const ptrdiff_t rs2 = Sizes<Mv::_rowsize,Mv::_colsize>::size;
        const ptrdiff_t rs = Sizes<rs1,rs2>::size;
        typedef typename Mu::cview_type Muv;
        typedef typename Vd::cview_type Vdv;
        typedef typename Ve::cview_type Vev;
        typedef typename Mv::cview_type Mvv;
        TMV_MAYBE_REF(Mu,Muv) Uv = U.cView();
        TMV_MAYBE_REF(Vd,Vdv) Dv = D.cView();
        TMV_MAYBE_REF(Ve,Vev) Ev = E.cView();
        TMV_MAYBE_REF(Mv,Mvv) Vv = V.cView();
        SVDecomposeFromBidiagonal_QR_Helper<-3,cs,rs,Muv,Vdv,Vev,Mvv>::call(
            Uv,Dv,Ev,Vv,UisI,VisI);
    }

    template <class Mu, class Vd, class Ve, class Mv>
    inline void SV_DecomposeFromBidiagonal_QR(
        BaseMatrix_Rec_Mutable<Mu>& U,
        BaseVector_Mutable<Vd>& D, BaseVector_Mutable<Ve>& E, 
        BaseMatrix_Rec_Mutable<Mv>& V, bool UisI, bool VisI)
    {
        TMVStaticAssert((Traits2<
                         typename Mu::value_type,
                         typename Vd::value_type>::samebase));
        TMVStaticAssert((Traits2<
                         typename Vd::value_type,
                         typename Ve::value_type>::sametype));
        TMVStaticAssert((Traits2<
                         typename Mu::value_type,
                         typename Mv::value_type>::samebase));
        TMVStaticAssert(Traits<typename Vd::value_type>::isreal);
        TMVStaticAssert(!Mu::_conj);
        TMVStaticAssert(!Mv::_conj);
        TMVStaticAssert((Sizes<Mu::_rowsize,Vd::_size>::same));
        TMVStaticAssert((Sizes<Vd::_size,IntTraits<Ve::_size>::Sp1>::same));
        TMVStaticAssert((Sizes<Mu::_rowsize,Mv::_colsize>::same));
        TMVStaticAssert((Sizes<Mu::_rowsize,Mv::_rowsize>::same));
        TMVAssert(D.size() == E.size()+1);
        if (U.cptr()) {
            TMVAssert(U.colsize() >= U.rowsize());
            TMVAssert(U.rowsize() == D.size());
        }
        if (V.cptr()) {
            TMVAssert(V.colsize() == D.size());
            TMVAssert(V.rowsize() >= V.colsize());
        }
        const ptrdiff_t cs = Mu::_colsize;
        const ptrdiff_t rs1 = Sizes<Mu::_rowsize,Vd::_size>::size;
        const ptrdiff_t rs2 = Sizes<Mv::_rowsize,Mv::_colsize>::size;
        const ptrdiff_t rs = Sizes<rs1,rs2>::size;
        typedef typename Mu::cview_type Muv;
        typedef typename Vd::cview_type Vdv;
        typedef typename Ve::cview_type Vev;
        typedef typename Mv::cview_type Mvv;
        TMV_MAYBE_REF(Mu,Muv) Uv = U.cView();
        TMV_MAYBE_REF(Vd,Vdv) Dv = D.cView();
        TMV_MAYBE_REF(Ve,Vev) Ev = E.cView();
        TMV_MAYBE_REF(Mv,Mvv) Vv = V.cView();
        SVDecomposeFromBidiagonal_QR_Helper<-2,cs,rs,Muv,Vdv,Vev,Mvv>::call(
            Uv,Dv,Ev,Vv,UisI,VisI);
    }

} // namespace tmv

#endif

