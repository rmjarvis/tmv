

#ifndef TMV_IntegerDet_H
#define TMV_IntegerDet_H

#include "TMV_BaseMatrix_Rec.h"

namespace tmv {

    // Defined in TMV_IntegerDet.cpp
    template <class T>
    T InstIntegerDet(const ConstMatrixView<T>& m);

    // This is a helper class that defines different operations for 
    // int and complex<int>.
    template <class T>
    struct Bareiss_Helper
    {
        typedef double double_type;
        typedef long double longdouble_type;

        static TMV_INLINE double_type to_double(T x)
        { return double_type(x); }

        static TMV_INLINE T convert(longdouble_type x) 
        { return x < 0. ? -T(floor(-x+0.5)) : T(floor(x+0.5)); }

        static TMV_INLINE bool isUnity(longdouble_type x) 
        { return TMV_ABS2(TMV_ABS2(x)-1.) < 0.1; }

        static TMV_INLINE bool isZero(longdouble_type x) 
        { return TMV_ABS2(x) < 0.1; }
    };

    template <class T>
    struct Bareiss_Helper<std::complex<T> >
    {
        typedef std::complex<T> CT;
        typedef std::complex<double> double_type;
        typedef std::complex<long double> longdouble_type;

        static TMV_INLINE double_type to_double(CT x)
        { return double_type(x.real(),x.imag()); }

        static TMV_INLINE CT convert(longdouble_type x) 
        {
            return CT(Bareiss_Helper<T>::convert(x.real()),
                      Bareiss_Helper<T>::convert(x.imag()));
        }

        static TMV_INLINE bool isUnity(longdouble_type x) 
        { return TMV_ABS2(TMV_ABS2(x) - 1.) < 0.1; }

        static TMV_INLINE bool isZero(longdouble_type x) 
        { return TMV_ABS2(x) < 0.1; }
    };

    template <int algo, class M>
    struct IntegerDet_Helper;

    // algo 11: Bareiss 1x1 algorithm
    template <class M>
    struct IntegerDet_Helper<11,M>
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
        typedef typename M::value_type T;
        static T call(const M& m)
        {
            typedef typename Bareiss_Helper<T>::longdouble_type DT;
            Matrix<DT,NoDivider> A(m);
            const int N = A.rowsize();

            // Do the 1x1 Bareiss step for each element down the diagonal.
            T det = 1; // Keep track of +-1 from row/col swaps.
            for (int k=0; k<N-1; ++k) {

                // Find the minimum non-zero value
                int imin = k, jmin = k;
                // Make sure we start with a non-zero value
                while (imin < N && 
                       Bareiss_Helper<T>::isZero(A.cref(imin,jmin))) ++imin;
                if (imin == N) return T(0);
                for (int j=k; j<N; ++j) {
                    if (Bareiss_Helper<T>::isUnity(A.cref(imin,jmin))) 
                        break;
                    for (int i=k; i<N; ++i) {
                        if (TMV_ABS2(A.cref(i,j))<TMV_ABS2(A.cref(imin,jmin)) &&
                            !Bareiss_Helper<T>::isZero(A.cref(i,j))) {
                            imin = i; jmin = j;
                            if (Bareiss_Helper<T>::isUnity(A.cref(imin,jmin))) 
                                break;
                        }
                    }
                }
                if (Bareiss_Helper<T>::isZero(A.cref(imin,jmin))) 
                    return T(0);

                if (imin != k) {
                    A.swapRows(imin,k);
                    det *= -1;
                }
                if (jmin != k) {
                    A.swapCols(jmin,k);
                    det *= -1;
                }

                // Now ready to do the 1x1 Bareiss step

                for (int j=k+1; j<N; ++j) for(int i=k+1; i<N; ++i) 
                    A.ref(i,j) = 
                        A.cref(k,k)*A.cref(i,j) -
                        A.cref(i,k)*A.cref(k,j);

                // If k > 0, then we have to divide by the previous step's
                // P value.
                if (k > 0) {
                    for (int j=k+1; j<N; ++j) for(int i=k+1; i<N; ++i) 
                        A.ref(i,j) /= A.cref(k-1,k-1);
                }
            }
            det *= Bareiss_Helper<T>::convert(A.cref(N-1,N-1));

            return det;
        }
    };

    // TODO: There is also a 2x2 (and higher nxn) Bareiss algorithm. 
    // Should check to see when this is faster.  

    // algo 90: Call inst
    template <class M>
    struct IntegerDet_Helper<90,M>
    {
        typedef typename M::value_type T;
        static TMV_INLINE T call(const M& m)
        { return InstIntegerDet(m.xView()); }
    };

    // algo 97: Conjugate
    template <class M>
    struct IntegerDet_Helper<97,M>
    {
        typedef typename M::value_type T;
        static TMV_INLINE T call(const M& m)
        {
            typedef typename M::const_conjugate_type Mc;
            Mc mc = m.conjugate();
            return TMV_CONJ(IntegerDet_Helper<-2,Mc>::call(mc));
        }
    };

    // algo -3: Determine which algorithm to use
    template <class M>
    struct IntegerDet_Helper<-3,M>
    {
        typedef typename M::value_type T;
        static TMV_INLINE T call(const M& m)
        {
            const int algo = 11;
            return IntegerDet_Helper<algo,M>::call(m);
        }
    };

    // algo -2: Check for inst
    template <class M>
    struct IntegerDet_Helper<-2,M>
    {
        typedef typename M::value_type T;
        static TMV_INLINE T call(const M& m)
        {
            const bool inst =
                (M::_colsize == TMV_UNKNOWN || M::_colsize > 16) &&
                (M::_rowsize == TMV_UNKNOWN || M::_rowsize > 16) &&
                Traits<T>::isinst;
            const int algo =
                M::_conj ? 97 :
                inst ? 90 :
                -3;
            return IntegerDet_Helper<algo,M>::call(m);
        }
    };

    template <class M>
    struct IntegerDet_Helper<-1,M>
    {
        typedef typename M::value_type T;
        static TMV_INLINE T call(const M& m)
        { return IntegerDet_Helper<-2,M>::call(m); }
    };

    template <class M>
    static inline typename M::value_type IntegerDet(
        const BaseMatrix_Rec<M>& m)
    {
        typedef typename M::value_type T;
        TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same));
        TMVAssert(m.colsize() == m.rowsize());
        TMVStaticAssert(Traits<T>::isinteger);
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return IntegerDet_Helper<-2,Mv>::call(mv);
    }

    template <class M>
    static inline typename M::value_type InlineIntegerDet(
        const BaseMatrix_Rec<M>& m)
    {
        typedef typename M::value_type T;
        TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same));
        TMVAssert(m.colsize() == m.rowsize());
        TMVStaticAssert(Traits<T>::isinteger);
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return IntegerDet_Helper<1,Mv>::call(mv);
    }

} // namespace tmv

#endif 
