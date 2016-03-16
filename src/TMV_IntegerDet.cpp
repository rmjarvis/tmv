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


#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_Permutation.h"
#include "tmv/TMV_PermutationArith.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_SimpleMatrix.h"
#include "TMV_IntegerDet.h"

#ifdef XDEBUG
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

    // This is a helper class that defines different operations for 
    // int and complex<int>.
    template <class T> 
    struct Helper
    {
        typedef double double_type;
        typedef long double longdouble_type;

        static double_type to_double(T x)
        { return double_type(x); }

        static T convert(longdouble_type x) 
        { return x < 0. ? -T(floor(-x+0.5)) : T(floor(x+0.5)); }

        static bool isUnity(longdouble_type x) 
        { return TMV_ABS2(TMV_ABS2(x)-1.) < 0.1; }

        static bool isZero(longdouble_type x) 
        { return TMV_ABS2(x) < 0.1; }
    };

    template <class T> 
    struct Helper<std::complex<T> >
    {
        typedef std::complex<T> CT;
        typedef std::complex<double> double_type;
        typedef std::complex<long double> longdouble_type;

        static double_type to_double(CT x)
        { return double_type(x.real(),x.imag()); }

        static CT convert(longdouble_type x) 
        { 
            return CT(Helper<T>::convert(x.real()),
                      Helper<T>::convert(x.imag()));
        }

        static bool isUnity(longdouble_type x) 
        { return TMV_ABS2(TMV_ABS2(x) - 1.) < 0.1; }

        static bool isZero(longdouble_type x) 
        { return TMV_ABS2(x) < 0.1; }
    };

    template <class T>
    static T Bareiss_1x1_IntegerDet(const GenMatrix<T>& A0)
    {
        // This algorithm is based on a paper by Erwin H. Bareiss:
        // Sylvester's Identity and Multistep Integer-Preserving
        // Gaussian Elimination, Mathematics of Computation,
        // Vol 22., No. 103, 1968, pp. 565-578.
        //
        // The basic idea of the algorithm is to partition the matrix A into
        // 
        // A = ( P  Q ) = ( P  0 ) ( I     P^-1 Q    )
        //     ( R  S )   ( R  I ) ( 0  S - R P^-1 Q )
        //
        // Then det(A) = det(P) det(S - R P^-1 Q).
        //
        // For this function, we'll take P to be a 1x1 submatrix.
        // (The next function will take up the case of P being 2x2, which 
        //  is faster, but this is still instructive.)
        // 
        // Then det(P) = P.
        //
        // And the P^-1 factor is a scalar, so it commutes with R and Q and 
        // we can factor it out:
        // det(S - R P^-1 Q) = P^-(n-1) det(P S - R Q)
        //
        // So, det(A) = P^-(n-2) det(P S - R Q)
        //
        // If n == 2, then we are done: det(A) = PS-RQ (as usual).
        // 
        // If n >= 3, then recursively calculate det(PS - RQ), which has
        // n reduced by 1.  Then divide by P^(n-2).
        //
        // Since we know det(A) has to be an exact integer, since it is the
        // determinant of a matrix of all integers, then we know that the 
        // det(PS-RQ) must be an exact multiple of P^(n-2), and integer
        // division is ok.
        //
        // The insight that Bareiss had was that we can actually do the 
        // division even sooner.  We can do the division at the next 
        // iteration of the algorithm when the submatrix is n-2 x n-2.
        // Then we can divide each element by P, which effectively divides
        // the determinant by P^(n-2).  
        //
        // This clearly work for real numbers.  Bareiss proved that it is ok
        // for integer division as well by proving that each element of 
        // that matrix is an exact multiple of P.
        //
        // I'll refer to the original paper for the full proof, but the 
        // basic idea is that each element in that matrix is P * the 
        // determinant of some 3x3 set of numbers from the original matrix A.
        // For example, the upper left element of that matrix has gone
        // through the same steps that would be performed to calculate just
        // the determinant of the upper left 3x3 submatrix of A.
        // Since that process ends with dividing by P, it is ok to divide it
        // by P at this point in the full determinant calculation.
        // Likewise each other element is the P * the 3x3 determinant of
        // ( A00 A01 A0j )
        // ( A10 A11 A1j )
        // ( Ai0 Ai1 Aij )
        // 
        // Supposedly, this process produces the smallest amount of growth
        // of the intermediate values.  Bareiss's proof of this seems 
        // somewhat weak, but there doesn't seem to be any refutation of 
        // this claim that I could find online, so I'll assume that it is 
        // true.
        //
        // However, the growth of intermediates is helped by three 
        // additional things we can do.
        //
        // First, at each step in the algorithm, we want P to be as small
        // as possible.  Even though we will be dividing it out a couple 
        // steps later, it helps to  have as small a number as possible for 
        // P in the P S calculation.
        // 
        // Second, we do all the calculations using long double storage.
        // This has 80 bits total, 64 of which are in the mantissa.
        // So it is quite a bit more precise than a 32-bit int, and has
        // the same precision as a 64-bit int.  However, in cases where
        // int would overflow in the intermediate values, long double will
        // often be able to salvage a useful answer.
        //
        // Also, this simplifies some of the math for complex numbers,
        // since division by complex<int> requires a special routine
        // to make sure the std library doesn't use abs to supposedly 
        // control rounding errors.
        //
         
        typedef typename Helper<T>::longdouble_type DT;
        SimpleMatrix<DT> A(A0);
        const ptrdiff_t N = A.rowsize();

        // Do the 1x1 Bareiss step for each element down the diagonal.
        T det = 1; // Keep track of +-1 from row/col swaps.
        for (ptrdiff_t k=0; k<N-1; ++k) {
            //cout<<"Start loop: k = "<<k<<endl;
            //cout<<"A = "<<A<<endl;

            // Find the minimum non-zero value
            ptrdiff_t imin = k, jmin = k;
            // Make sure we start with a non-zero value
            while (imin < N && Helper<T>::isZero(A.cref(imin,jmin))) ++imin;
            if (imin == N) return T(0);
            for (ptrdiff_t j=k; j<N; ++j) {
                if (Helper<T>::isUnity(A.cref(imin,jmin)))  break;
                for (ptrdiff_t i=k; i<N; ++i) {
                    if (TMV_ABS2(A.cref(i,j)) < TMV_ABS2(A.cref(imin,jmin)) &&
                        !Helper<T>::isZero(A.cref(i,j))) {
                        imin = i; jmin = j;
                        if (Helper<T>::isUnity(A.cref(imin,jmin)))  break;
                    }
                }
            }
            //cout<<"min value is "<<A.cref(imin,jmin)<<" at "<<imin<<','<<jmin<<endl;
            if (Helper<T>::isZero(A.cref(imin,jmin))) return T(0);

            if (imin != k) {
                A.swapRows(imin,k);
                det *= -1;
            }
            if (jmin != k) {
                A.swapCols(jmin,k);
                det *= -1;
            }

            // Now ready to do the 1x1 Bareiss step
            //cout<<"P = "<<A.cref(k,k)<<endl;

            for (ptrdiff_t j=k+1; j<N; ++j) for(ptrdiff_t i=k+1; i<N; ++i) 
                A.ref(i,j) = A.cref(k,k)*A.cref(i,j) - A.cref(i,k)*A.cref(k,j);

            // If k > 0, then we have to divide by the previous step's
            // P value.
            if (k > 0) {
                for (ptrdiff_t j=k+1; j<N; ++j) for(ptrdiff_t i=k+1; i<N; ++i) 
                    A.ref(i,j) /= A.cref(k-1,k-1);
            }
        }
        det *= Helper<T>::convert(A.cref(N-1,N-1));
        //cout<<"Final det = "<<det<<endl;

        return det;
    }

    template <class T>
    T IntegerDet(const GenMatrix<T>& A)
    {
#ifdef XDEBUG
        typedef typename Helper<T>::double_type DT;
        Matrix<DT> D(A);
        DT det_act = D.det();
        cout<<"Start IntegerDet: A = "<<A<<endl;
        cout<<"det_act = "<<det_act<<endl;
#endif
        TMVAssert(A.isSquare());

        const ptrdiff_t N = A.colsize();
        T det = 
            N == 0 ? T(1) :
            N == 1 ? A.cref(0,0) :
            N == 2 ? (A.cref(0,0) * A.cref(1,1) - A.cref(0,1) * A.cref(1,0)) :
            N == 3 ? (
                A.cref(0,0)*(A.cref(1,1)*A.cref(2,2)-A.cref(1,2)*A.cref(2,1)) -
                A.cref(0,1)*(A.cref(1,0)*A.cref(2,2)-A.cref(1,2)*A.cref(2,0)) +
                A.cref(0,2)*(A.cref(1,0)*A.cref(2,1)-A.cref(1,1)*A.cref(2,0)) 
            ) :
            Bareiss_1x1_IntegerDet(A);

#ifdef XDEBUG
        cout<<"Done: det = "<<det<<endl;
        DT det_d = Helper<T>::to_double(det);
        cout<<"abs(diff) = "<<std::abs(det_d-det_act)<<endl;
        if (std::abs(det_d - det_act) > 1.e-4) {
            cerr<<"IntegerDet: A = "<<A<<endl;
            cerr<<"det_act = "<<det_act<<endl;
            cerr<<"det = "<<det<<endl;
            cerr<<"det_d = "<<det_d<<endl;
            abort();
        }
#endif
        return det;
    }

#define InstFile "TMV_IntegerDet.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


