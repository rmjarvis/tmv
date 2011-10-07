

//
// This file contains the code to implement Householder reflections.
//
// A Householder reflection is a unitary matrix H which transforms
// a given vector x into y e0 = (y,0,0,0...0) where |y| = |x|
//
// H is of the form: 
//
// H = I - beta u ut = I - beta (x - y e0) (xt - y* e0t)
//
// H x = (1 - beta (|x|^2 - y* x0)) x + beta (y |x|^2 - |y|^2 x0) e0
//
// Since H x is supposed to equal y e0, the coefficient of x must be 0, so
//
// beta = 1 / (|x|^2 - y* x0)
//
// With this definition for beta, the coefficient of e0 is:
//
// (y |x|^2 - |y|^2 x0) / (|x|^2 - y* x0)
// = |y|^2 (y - x0) / y* (y - x0) = y
//
// Check that H is unitary:
//
// HtH = (I - beta* (x - y e0) (xt - y* e0t))(I - beta (x - y e0) (xt - y* e0t))
//     = I - (beta*+beta) (x xt - y e0 xt - y* x e0t + |y|^2 e0 e0t)
//         + |beta|^2 (x - y e0) (|x|^2 - y* x0 - y x0* + |y|^2) (xt - y* e0t)
//     = I - [beta + beta* - |beta|^2 (1/beta + 1/beta*)] (x-ye0)(xt-y* e0t)
//     = I - (beta + beta* - beta* - beta)(x-y e0)(xt-y* e0t)
//     = I
//
// Note that we still have a choice on the value for y.  The only constraint
// is that |y| = |x| (the vector norm).  So y = |x| exp(i theta)
//
// There are two options which make some sense: 
// choosing a real y ( +-|x| ) or choosing y so that beta is real
// ( 1/(|x|^2 +- |x||x0|) ).
//
// The better choice is to make beta real, since then the matrix is
// Hermitian as well as unitary (ie. H = Ht, H*H = I), which also
// leads to the determinant being -1, which is actually required for
// it to be a "reflection" matrix.  It also has the advantage of (slightly)
// quicker multiplies, since the multiplication of H involves a multiplication
// by beta.  This step takes half the time if beta is real.  Of course, it's
// not the bulk of the calculation for large matrices, hence the "(slightly)".
//
// Unfortunately, LaPack made the wrong choice here.  They chose to make
// y real.  The _only_ advantage (that I can see) for this choice is 
// in calculating the SVD of a matrix.  Then the bidiagonalization
// using Householder matrices goes directly to a real bidiagonal, so
// it doesn't require a second step using Givens matrices to make it real.
// This is a pretty minor advantage which doesn't offset the speed loss
// from the complex multiplies for the betas, and makes keeping track of 
// the determinant much harder.
//
// Prior to TMV v0.70, I chose the LaPack convention to keep compatibility
// with their routines.  However, with TMV v0.70, I decided to switch to 
// using the true reflection convention, so now I have real beta values.
//
// We still have two possible real values for beta:
// ( 1/(|x|^2 +- |x||x0|) ).
//
// We want to minimize the effects of rounding errors, which means maximize
// 1/beta = |x|^2 +- |x| |x0|
// Clearly the + option is larger, so that is the choice we use.
//
// The corresponding value of y is then: y = -|x| (x0/|x0|) 
// (If x0 == 0, use y = |x|.)
//
//
// The determinant of this Householder non-reflection is:
//
// det(H) = det( I - beta (x+|x|e0) (x+|x|e0)t )
//
// A simple, but not (to me) obvious identity is:
//
// det(I-alpha u ut) = 1 - alpha |u|^2
//
// Aside: The proof is pretty cool and comes from the following identity:
// (the so-called Matrix Determinant Lemma)
//   ( I   0 ) ( I+uvT  u ) (  I   0 ) = ( I    u   )
//   ( vT  1 ) (   0    1 ) ( -vT  1 )   ( 0  1+vTu )
// The determinant of the left side is det(I+uvT).
// The determinant of the right side is 1+vTu.
// Our special case comes from v = -alpha u.conj().
//
//
// Thus, det(H) = 1 - beta |x - y e0|^2 
// |x - y e0|^2 = |x0-y|^2 + |x_1..N|^2
//              = |x0|^2 + |x_1..N|^2 - x0* y - x0 y* + |y|^2
//              = |x|^2 - x0* y - x0 y* + |x|^2
//              = 2|x|^2 - 2 Re( x0* (-|x| x0/|x0|) )
//              = 2|x|^2 + 2 Re( |x0| |x| )
//              = 2|x|^2 + 2|x0||x|
// det(H) = 1 - 2 beta (|x|^2 + |x0||x|)
//        = 1 - 2 (|x|^2 + |x0||x|) / (|x|^2 + |x0||x|)
//        = -1
//
// There are two further complications to deal with in practice.
//
// First, we don't introduce new storage to store the u part of the 
// Householder matrix definition.  Instead, we store it in the same location
// as the original vector that we reflect into y e0.  However, this means
// that we don't have space to store the first element of u, since we
// need to keep the resulting y value there.  So the solution is to 
// recognize that there is another degeneracy that we can take advantage of.
// Specifically, if u is scaled by a value w, and beta is scaled by 1/w^2,
// then the Householder matrix is unchanged.
//
// This means that we can scale u so that the first element is some given
// value.  This is valid, since u0 can never be zero:  
// If x0 != 0:
//     u0 = x0 - y = x0 + |x| x0 / |x0| = x0 * (1 + |x|/|x0|) != 0
// If x0 == 0:
//     u0 = -y = -|x| != 0  (assuming |x| != 0, see below)
//
// So we can scale u by 1/u0, and get a version that has u0 == 1.
// In addition we also have to scale beta by u0^2.  
// Then all further usage of the Householder matrix can proceed
// with the implicit value that u(0) == 1, rather than the value that
// is actually stored in that location.
//
// Second, the above prescription runs into a problem in the case of an
// exactly zero input vector, since x0 == 0 and |x| == 0, so u0 = 0.
// Furthermore, if x0 != 0, but the rest of x is zero, then there is 
// no point in doing the reflection, since the vector is already in the
// right form.  So we can save some calculation by using H = I instead
// of the normal calculation.  (In this case, H is not a reflection, so 
// det = 1, not -1.)  We mark this case by setting beta = 0.
//

#ifndef TMV_Householder_H
#define TMV_Householder_H

#include "tmv/TMV_BaseVector.h"
#include "tmv/TMV_BaseMatrix_Rec.h"
#include "tmv/TMV_BaseMatrix_Tri.h"
#include "tmv/TMV_MultMV.h"
#include "tmv/TMV_Rank1VVM.h"

#ifdef PRINTALGO_HOUSE
#include <iostream>
#include "TMV_MatrixIO.h"
#include "TMV_VectorIO.h"
#endif

//#define XDEBUG_HOUSE

#ifdef XDEBUG_HOUSE
#include <iostream>
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_AddMM.h"
#include "tmv/TMV_SumMX.h"
#include "tmv/TMV_NormM.h"
#include "tmv/TMV_NormV.h"
#include "tmv/TMV_Norm.h"
#include "tmv/TMV_MatrixIO.h"
#include "tmv/TMV_VectorIO.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_CopyM.h"
using std::cerr;
using std::endl;
#endif

namespace tmv {

    // Defined in TMV_Householder.cpp
    template <class T>
    void InstHouseholderReflect(
        T& x0, VectorView<T> u, typename Traits<T>::real_type& beta);
    template <class T>
    bool InstHouseholderUnReflect(
        T& x0, VectorView<T> u, typename Traits<T>::real_type& beta);
    template <class T>
    void InstHouseholderUnpack(
        VectorView<T> u, typename Traits<T>::real_type& beta, T& u0);
    template <class T1, int C1, class T2>
    void InstHouseholderMultEq(
        const ConstVectorView<T1,C1>& u, typename Traits<T1>::real_type beta,
        VectorView<T2> m0, MatrixView<T2> mx, VectorView<T2> temp);
    template <class T>
    void InstBlockHouseholderAugment(
        const ConstMatrixView<T>& Y, UpperTriMatrixView<T> Z,
        typename Traits<T>::real_type beta);
    template <class T>
    void InstBlockHouseholderMakeZ(
        const ConstMatrixView<T>& Y, UpperTriMatrixView<T> Z,
        const ConstVectorView<typename Traits<T>::real_type>& beta);
    template <class T1, int C1, class T2>
    void InstBlockHouseholderLMult(
        const ConstMatrixView<T1,C1>& Y,
        const ConstUpperTriMatrixView<T1,C1>& Z,
        MatrixView<T2> m2, MatrixView<T2> temp);
    template <class T1, int C1, class T2>
    void InstBlockHouseholderLDiv(
        const ConstMatrixView<T1,C1>& Y,
        const ConstUpperTriMatrixView<T1,C1>& Z,
        MatrixView<T2> m2, MatrixView<T2> temp);
    template <class T>
    void InstBlockHouseholderUnpack(
        MatrixView<T> Y, const ConstUpperTriMatrixView<T>& Z,
        UpperTriMatrixView<T> temp);
    template <class T>
    void InstBlock2HouseholderAugment(
        const ConstMatrixView<T>& Y, UpperTriMatrixView<T> Z,
        typename Traits<T>::real_type beta);
    template <class T>
    void InstBlock2HouseholderMakeZ(
        const ConstMatrixView<T>& Y, UpperTriMatrixView<T> Z,
        const ConstVectorView<typename Traits<T>::real_type>& beta);
    template <class T1, int C1, class T2, int C2>
    void InstBlock2HouseholderLMult(
        const ConstMatrixView<T1,C1>& Y,
        const ConstUpperTriMatrixView<T1,C1>& Z,
        MatrixView<T2,C2> ma, MatrixView<T2> mb, MatrixView<T2> temp);
    template <class T1, int C1, class T2, int C2>
    void InstBlock2HouseholderLDiv(
        const ConstMatrixView<T1,C1>& Y,
        const ConstUpperTriMatrixView<T1,C1>& Z,
        MatrixView<T2,C2> ma, MatrixView<T2> mb, MatrixView<T2> temp);



    //
    // HouseholderReflect
    // Finds the Householder matrix H which rotates v into y e0.
    // The vector v of the Householder matrix is stored in v,
    // except for the first element.
    // Beta is the return value.
    //
    
    template <int algo, class V>
    struct HouseholderReflect_Helper;

    // algo 11: Normal case
    template <class V>
    struct HouseholderReflect_Helper<11,V>
    {
        typedef typename V::value_type T;
        typedef typename V::real_type RT;
        static void call(T& x0, V& u, RT& beta)
        {
#ifdef XDEBUG_HOUSE
            Vector<T> xx(u.size()+1);
            xx(0) = x0;
            xx.subVector(1,xx.size()) = u;
#endif
            // Check for possibility of an early exit.
            if (u.size() == 0) {
                beta = RT(0);
                return;
            }

            // Since this routine involves squares of elements, we risk overflow
            // and underflow problems if done naively.  The simplest (although
            // probably not the most efficient) solution is to scale all the 
            // intermediate values by the maximum abs value in the vector.

            RT scale = u.maxAbs2Element();
            RT absx0 = TMV_ABS(x0);
            if (absx0 > scale) scale = absx0;

            if (TMV_Underflow(scale)) {
                // Then the situation is hopeless, and we should just zero out
                // the whole vector.
                u.setZero();
                x0 = T(0);
                beta = RT(0);
                return;
            }

            // Determine normx = |x|
            RT invscale = RT(1)/scale;
            RT normsqx = u.normSq(invscale);

            // if all of x other than first element are 0, H is identity
            if (normsqx == RT(0)) {
                // Set them all to x explicitly in case underflow led to the 0.
                u.setZero();
                beta = RT(0);
                return;
            }

            // Finish calculating normx in usual case
            absx0 *= invscale;
            x0 *= invscale;
            RT normsqx0 = absx0*absx0;
            normsqx += normsqx0;
            RT normx = TMV_SQRT(normsqx);

            // y = -|x| (x0/|x0|) 
            T y = absx0 == RT(0) ? normx : -normx*x0/absx0;

            // beta = 1 / (|x|^2 + |x||x0|)
            // H = I - beta v vt
            // with u = x - y e0 in first column
            // Renormalize beta,u so that u(0) = 1
            T u0 = x0-y;
            RT normu0 = TMV_NORM(u0);
            beta = normu0 / (normsqx + normx * absx0);
            T invu0 = RT(1)/u0;

            // Sometimes this combination can underflow, so check.
            T scaled_invu0 = invu0 * invscale;
            if ((TMV_REAL(invu0)!=RT(0) && TMV_Underflow(TMV_REAL(scaled_invu0))) ||
                (TMV_IMAG(invu0)!=RT(0) && TMV_Underflow(TMV_IMAG(scaled_invu0)))) {
                u *= invscale;
                u *= invu0;
            } else {
                u *= scaled_invu0;
            }

            x0 = y*scale;

#ifdef XDEBUG_HOUSE
            Vector<T> uu(xx.size());
            uu(0) = T(1);
            uu.subVector(1,uu.size()) = u;
            Matrix<T> H = T(1)-beta*(uu^uu.conjugate());
            // Check the following:
            // H * xx = (y,0,0.0...)
            // |y| = Norm(xx)
            Vector<T> Hxx = H * xx;
            if (TMV_ABS(TMV_ABS(Hxx(0))-TMV_ABS(x0)) > 0.0001*TMV_ABS(x0) ||
                TMV_ABS(Norm(xx)-TMV_ABS(x0)) > 0.0001*TMV_ABS(x0) ||
                Norm(Hxx.subVector(1,Hxx.size())) > 0.0001*TMV_ABS(x0)) {
                cerr<<"Householder Reflect:\n";
                cerr<<"Input: x = "<<xx<<endl;
                cerr<<"x0 = "<<x0<<endl;
                cerr<<"Norm(x) = "<<Norm(xx)<<endl;
                cerr<<"Output: uu= "<<uu<<endl;
                cerr<<"beta = "<<beta<<endl;
                cerr<<"H = "<<H<<endl;
                cerr<<"xx = "<<xx<<endl;
                cerr<<"Hxx = "<<Hxx<<endl;
                cerr<<"abs(abs(hxx(0))-abs(x0)) = "<<
                    TMV_ABS(TMV_ABS(Hxx(0))-TMV_ABS(x0))<<endl;
                cerr<<"abs(Norm(xx)-abs(x0)) = "<<
                    TMV_ABS(Norm(xx))-TMV_ABS(x0)<<endl;
                cerr<<"Norm(Hxx(1,N)) = "<<Norm(Hxx.subVector(1,Hxx.size()))<<endl;
                abort();
            }
#endif
        }
    };

    // algo 90: Call inst
    template <class V>
    struct HouseholderReflect_Helper<90,V>
    {
        typedef typename V::value_type T;
        typedef typename V::real_type RT;
        static void call(T& x0, V& u, RT& beta)
        { InstHouseholderReflect(x0,u.xView(),beta); }
    };

    // algo 97: Conjugate u
    template <class V>
    struct HouseholderReflect_Helper<97,V>
    {
        typedef typename V::value_type T;
        typedef typename V::real_type RT;
        static void call(T& x0, V& u, RT& beta)
        { 
            typedef typename V::conjugate_type Vc;
            Vc uc = u.conjugate();
            x0 = TMV_CONJ(x0);
            HouseholderReflect_Helper<-2,Vc>::call(x0,uc,beta); 
            x0 = TMV_CONJ(x0);
        }
    };

    // algo -3: Select algorithm
    template <class V>
    struct HouseholderReflect_Helper<-3,V>
    {
        typedef typename V::value_type T;
        typedef typename V::real_type RT;
        static void call(T& x0, V& u, RT& beta)
        { 
            const int algo = 11;
#ifdef PRINTALGO_HOUSE
            std::cout<<"Inline HouseholderReflect\n";
            std::cout<<"x0 = "<<x0<<std::endl;
            std::cout<<"u = "<<u<<std::endl;
            std::cout<<"beta = "<<beta<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            HouseholderReflect_Helper<algo,V>::call(x0,u,beta); 
#ifdef PRINTALGO_HOUSE
            std::cout<<"x0 => "<<x0<<std::endl;
            std::cout<<"u => "<<u<<std::endl;
            std::cout<<"beta => "<<beta<<std::endl;
#endif
        }
    };

    // algo -2: Check for inst
    template <class V>
    struct HouseholderReflect_Helper<-2,V>
    {
        typedef typename V::value_type T;
        typedef typename V::real_type RT;
        static void call(T& x0, V& u, RT& beta)
        { 
            const bool inst = 
                (V::_size == TMV_UNKNOWN || V::_size > 16) &&
                Traits<T>::isinst;
            const int algo = 
                V::_conj ? 97 : 
                inst ? 90 :
                -3;
            HouseholderReflect_Helper<algo,V>::call(x0,u,beta); 
        }
    };

    template <class V>
    struct HouseholderReflect_Helper<-1,V>
    {
        typedef typename V::value_type T;
        typedef typename V::real_type RT;
        static void call(T& x0, V& u, RT& beta)
        { HouseholderReflect_Helper<-2,V>::call(x0,u,beta); }
    };

    template <class V>
    static inline void InlineHouseholderReflect(
        typename V::value_type& x0, BaseVector_Mutable<V>& u, 
        typename V::real_type& beta)
    {
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(V,Vv) uv = u.cView();
        HouseholderReflect_Helper<-3,Vv>::call(x0,uv,beta);
    }

    template <class V>
    static inline void HouseholderReflect(
        typename V::value_type& x0, BaseVector_Mutable<V>& u, 
        typename V::real_type& beta)
    {
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(V,Vv) uv = u.cView();
        HouseholderReflect_Helper<-2,Vv>::call(x0,uv,beta);
    }

    // x0 is also allowed to be a ConjRef
    template <class V>
    static inline void HouseholderReflect(
        ConjRef<typename V::value_type> x0, BaseVector_Mutable<V>& u, 
        typename V::real_type& beta)
    {
        typename V::value_type temp = x0;
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(V,Vv) uv = u.cView();
        HouseholderReflect_Helper<-2,Vv>::call(temp,uv,beta);
        x0 = temp;
    }


    //
    // HouseholderUnReflect
    // This is similar, except that the roles of y and x0 are swapped.
    // That is, the final rotated vector y e0 is presumed known, as is 
    // the bulk of the input vector.  
    // The unreflected value of the first element in the vector is returned 
    // as y.
    // The rest of the vector is transformed into the Householder vector.
    // This is used for downdating QR decompositions.
    // The return value is true if successful, false if |y| < |x|.
    //
 
    template <int algo, class V>
    struct HouseholderUnReflect_Helper;

    // algo 11: Normal case
    template <class V>
    struct HouseholderUnReflect_Helper<11,V>
    {
        typedef typename V::value_type T;
        typedef typename V::real_type RT;
        static bool call(T& y, V& u, RT& beta)
        {
            typedef typename V::value_type T;
            typedef typename V::real_type RT;

            RT scale = u.maxAbs2Element();
            RT absy = TMV_ABS(y);
            if (absy > scale) scale = absy;

            if (TMV_Underflow(scale)) {
                // Then the situation is hopeless, and we should just zero out
                // the whole vector.
                u.setZero();
                y = T(0);
                beta = RT(0);
                return true; // But still declared successful.
            }

            // Determine |y| = |x| now, but need to find x0 from |x1| and |y|.
            RT invscale = RT(1)/scale;
            RT normsqx1 = u.normSq(invscale);

            // if all of x other than first element are 0, H is identity
            if (normsqx1 == RT(0)) {
                // Set them all to x explicitly in case underflow led to the 0.
                u.setZero();
                beta = RT(0);
                return true;
            }

            // Finish calculating x0 in usual case
            absy *= invscale;
            y *= invscale;
            RT normsqx = absy*absy;
            RT normsqx0 = normsqx - normsqx1;
            if (normsqx0 < RT(0)) return false;

            // y = -|x| (x0/|x0|) 
            // x0 = -|x0| (y/|y|)
            RT absx0 = TMV_SQRT(normsqx0);
            T x0 = -absx0*y/absy;

            // beta = 1 / (|x|^2 + |x||x0|)
            // H = I - beta v vt
            // with u = x - y e0 in first column
            // Renormalize beta,u so that u(0) = 1
            T u0 = x0-y;
            RT normu0 = TMV_NORM(u0);
            beta = normu0 / (normsqx + absy * absx0);
            T invu0 = RT(1)/u0;

            // Sometimes this combination can underflow, so check.
            T scaled_invu0 = invu0 * invscale;
            if ((TMV_REAL(invu0)!=RT(0) && TMV_Underflow(TMV_REAL(scaled_invu0))) ||
                (TMV_IMAG(invu0)!=RT(0) && TMV_Underflow(TMV_IMAG(scaled_invu0)))) {
                u *= invscale;
                u *= invu0;
            } else {
                u *= scaled_invu0;
            }

            y = x0*scale;
            return true;
        }
    };

    // algo 90: Call inst
    template <class V>
    struct HouseholderUnReflect_Helper<90,V>
    {
        typedef typename V::value_type T;
        typedef typename V::real_type RT;
        static bool call(T& y, V& u, RT& beta)
        { return InstHouseholderUnReflect(y,u.xView(),beta); }
    };

    // algo 97: Conjugate u
    template <class V>
    struct HouseholderUnReflect_Helper<97,V>
    {
        typedef typename V::value_type T;
        typedef typename V::real_type RT;
        static bool call(T& y, V& u, RT& beta)
        { 
            typedef typename V::conjugate_type Vc;
            Vc uc = u.conjugate();
            y = TMV_CONJ(y);
            return HouseholderUnReflect_Helper<-2,Vc>::call(y,uc,beta); 
            y = TMV_CONJ(y);
        }
    };

    // algo -3: Select algorithm
    template <class V>
    struct HouseholderUnReflect_Helper<-3,V>
    {
        typedef typename V::value_type T;
        typedef typename V::real_type RT;
        static bool call(T& y, V& u, RT& beta)
        { 
            const int algo = 11;
#ifdef PRINTALGO_HOUSE
            std::cout<<"Inline HouseholderUnReflect\n";
            std::cout<<"y = "<<y<<std::endl;
            std::cout<<"u = "<<u<<std::endl;
            std::cout<<"beta = "<<beta<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            return HouseholderUnReflect_Helper<algo,V>::call(y,u,beta); 
#ifdef PRINTALGO_HOUSE
            std::cout<<"y => "<<y<<std::endl;
            std::cout<<"u => "<<u<<std::endl;
            std::cout<<"beta => "<<beta<<std::endl;
#endif
        }
    };

    // algo -2: Check for inst
    template <class V>
    struct HouseholderUnReflect_Helper<-2,V>
    {
        typedef typename V::value_type T;
        typedef typename V::real_type RT;
        static bool call(T& y, V& u, RT& beta)
        { 
            const bool inst = 
                (V::_size == TMV_UNKNOWN || V::_size > 16) &&
                Traits<T>::isinst;
            const int algo = 
                V::_conj ? 97 : 
                inst ? 90 :
                -3;
            return HouseholderUnReflect_Helper<algo,V>::call(y,u,beta); 
        }
    };

    template <class V>
    struct HouseholderUnReflect_Helper<-1,V>
    {
        typedef typename V::value_type T;
        typedef typename V::real_type RT;
        static bool call(T& y, V& u, RT& beta)
        { return HouseholderUnReflect_Helper<-2,V>::call(y,u,beta); }
    };

    template <class V>
    static inline bool InlineHouseholderUnReflect(
        typename V::value_type& y, BaseVector_Mutable<V>& u, 
        typename V::real_type& beta)
    {
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(V,Vv) uv = u.cView();
        return HouseholderUnReflect_Helper<-3,Vv>::call(y,uv,beta);
    }

    template <class V>
    static inline bool HouseholderUnReflect(
        typename V::value_type& y, BaseVector_Mutable<V>& u, 
        typename V::real_type& beta)
    {
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(V,Vv) uv = u.cView();
        return HouseholderUnReflect_Helper<-2,Vv>::call(y,uv,beta);
    }

    // y is also allowed to be a ConjRef
    template <class V>
    static inline bool HouseholderUnReflect(
        ConjRef<typename V::value_type> y, BaseVector_Mutable<V>& u, 
        typename V::real_type& beta)
    {
        typename V::value_type temp = y;
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(V,Vv) uv = u.cView();
        bool ret = HouseholderUnReflect_Helper<-2,Vv>::call(temp,uv,beta);
        y = temp;
        return ret;
    }


    //
    // HouseholderUnpack
    // This routine takes u <- H*e0
    //
    
    template <int algo, class V>
    struct HouseholderUnpack_Helper;
    
    // algo 11: Normal case
    template <class V>
    struct HouseholderUnpack_Helper<11,V>
    {
        typedef typename V::reference Ref;
        typedef typename V::value_type T;
        typedef typename V::real_type RT;
        static void call(V& u, RT beta, Ref u0)
        {
            if (beta == T(0)) {
                // then all but u(0) is already 0
                u0 = T(1);
            } else {
                // u <- (I-beta (1 u) (1 u)t) e0 = e0 - beta (1 u)
                u0 = T(1)-beta;
                u *= -beta;
            }
        }
    };

    // algo 90: Call inst
    template <class V>
    struct HouseholderUnpack_Helper<90,V>
    {
        typedef typename V::reference Ref;
        typedef typename V::value_type T;
        typedef typename V::real_type RT;
        static void call(V& u, RT beta, Ref u0)
        { InstHouseholderUnpack(u.xView(),beta,u0); }
    };

    // algo 97: Conjugate u
    template <class V>
    struct HouseholderUnpack_Helper<97,V>
    {
        typedef typename V::reference Ref;
        typedef typename V::value_type T;
        typedef typename V::real_type RT;
        static void call(V& u, RT beta, Ref u0)
        { 
            typedef typename V::conjugate_type Vc;
            Vc uc = u.conjugate();
            HouseholderUnpack_Helper<-2,Vc>::call(uc,beta,TMV_CONJ(u0)); 
        }
    };

    // algo -3: Select algorithm
    template <class V>
    struct HouseholderUnpack_Helper<-3,V>
    {
        typedef typename V::reference Ref;
        typedef typename V::value_type T;
        typedef typename V::real_type RT;
        static void call(V& u, RT beta, Ref u0)
        { 
            const int algo = 11;
#ifdef PRINTALGO_HOUSE
            std::cout<<"Inline HouseholderUnPack\n";
            std::cout<<"u = "<<u<<std::endl;
            std::cout<<"beta = "<<beta<<std::endl;
            std::cout<<"u0 = "<<u0<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            HouseholderUnpack_Helper<algo,V>::call(u,beta,u0); 
#ifdef PRINTALGO_HOUSE
            std::cout<<"u => "<<u<<std::endl;
            std::cout<<"u0 = >"<<u0<<std::endl;
#endif
        }
    };

    // algo -2: Check for inst
    template <class V>
    struct HouseholderUnpack_Helper<-2,V>
    {
        typedef typename V::reference Ref;
        typedef typename V::value_type T;
        typedef typename V::real_type RT;
        static void call(V& u, RT beta, Ref u0)
        { 
            const bool inst = 
                (V::_size == TMV_UNKNOWN || V::_size > 16) &&
                Traits<T>::isinst;
            const int algo = 
                V::_conj ? 97 : 
                inst ? 90 :
                -3;
            HouseholderUnpack_Helper<algo,V>::call(u,beta,u0); 
        }
    };

    template <class V>
    struct HouseholderUnpack_Helper<-1,V>
    {
        typedef typename V::reference Ref;
        typedef typename V::value_type T;
        typedef typename V::real_type RT;
        static void call(V& u, RT beta, Ref u0)
        { HouseholderUnpack_Helper<-2,V>::call(u,beta,u0); }
    };

    template <class V>
    static inline void InlineHouseholderUnpack(
        BaseVector_Mutable<V>& u, 
        typename V::real_type beta, typename V::reference u0)
    {
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(V,Vv) uv = u.cView();
        HouseholderUnpack_Helper<-3,Vv>::call(uv,beta,u0);
    }

    template <class V>
    static inline void HouseholderUnpack(
        BaseVector_Mutable<V>& u, 
        typename V::real_type beta, typename V::reference u0)
    {
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(V,Vv) uv = u.cView();
        HouseholderUnpack_Helper<-2,Vv>::call(uv,beta,u0);
    }



    //
    // HouseholderMultEq
    // This routine takes 
    // ( m0 ) <- H ( m0 )
    // ( mx )      ( mx )
    //

    template <int algo, class V, class M0, class Mx, class Vt>
    struct HouseholderMultEq_Helper;

    // algo 11: Normal case
    template <class V, class M0, class Mx, class Vt>
    struct HouseholderMultEq_Helper<11,V,M0,Mx,Vt>
    {
        typedef typename V::real_type RT;
        static void call(const V& u, RT beta, M0& m0, Mx& mx, Vt& temp)
        {
            // ( m0 ) <- H ( m0 ) = [ ( m0 ) - beta ( 1 ) ( 1 ut ) ( m0 ) ]
            // ( mx )      ( mx )   [ ( mx )        ( u )          ( mx ) ]
            // 
            // ( m0 ) -= beta (   m0 + ut mx   )
            // ( mx )         ( u (m0 + ut mx) )

            TMVAssert(u.size() == mx.colsize());
            TMVAssert(m0.size() == mx.rowsize());
            TMVAssert(temp.size() == mx.rowsize());
#ifdef XDEBUG_HOUSE
            std::cout<<"Start Householder::multEq"<<std::endl;
            std::cout<<"u = "<<u<<std::endl;
            std::cout<<"beta = "<<beta<<std::endl;
            std::cout<<"m0 = "<<m0.vec()<<std::endl;
            std::cout<<"mx = "<<mx.mat()<<std::endl;
            typedef typename V::value_type T;
            typedef typename Mx::value_type T2;
            Matrix<T2> mm(mx.colsize()+1,m0.size());
            mm.row(0) = m0.vec();
            mm.rowRange(1,mm.colsize()) = mx.mat();
#endif
            typedef typename V::real_type RT;

            if (m0.size() > 0 && beta != RT(0)) {
                //std::cout<<"Not trivial.\n";
                temp = u.conjugate() * mx.mat();
                //std::cout<<"temp = "<<temp<<std::endl;
                temp += m0.vec();
                //std::cout<<"temp => "<<temp<<std::endl;
                temp *= -beta;
                //std::cout<<"temp => "<<temp<<std::endl;

                m0.vec() += temp;
                //std::cout<<"m0 => "<<m0<<std::endl;
                mx.mat() += u^temp;
                //std::cout<<"m0 => "<<mx<<std::endl;
            }
#ifdef XDEBUG_HOUSE
            Vector<T> uu(u.size()+1);
            uu(0) = T(1);
            uu.subVector(1,uu.size()) = u;
            Matrix<T> H = T(1) - beta*(uu^uu.conjugate());
            Matrix<T2> Hm = H * mm;
            Matrix<T2> Hm2(mx.colsize()+1,m0.size());
            Hm2.row(0) = m0.vec();
            Hm2.rowRange(1,Hm2.colsize()) = mx.mat();
            if (!(Norm(Hm-Hm2) <= 0.001*Norm(Hm))) {
                cerr<<"Householder::multEq\n";
                cerr<<"Input: m = "<<mm<<endl;
                cerr<<"uu = "<<uu<<endl;
                cerr<<"beta = "<<beta<<endl;
                cerr<<"H = "<<H<<endl;
                cerr<<"Hm = "<<Hm<<endl;
                cerr<<"Output: m = "<<Hm2<<endl;
                abort();
            }
            //std::cout<<"mx => "<<mx.mat()<<std::endl;
#endif
        }
    };

    // algo 82: Use temp for m0
    // This would only be necessary if m0 and mx come from different
    // matrices and one of them is conjugated.  Since usually m0 and mx
    // are part of the same actual matrix, this wouldn't usually get called.
    template <class V, class M0, class Mx, class Vt>
    struct HouseholderMultEq_Helper<82,V,M0,Mx,Vt>
    {
        typedef typename V::real_type RT;
        static void call(const V& u, RT beta, M0& m0, Mx& mx, Vt& temp)
        { 
            typedef typename M0::copy_type M0c;
            typedef typename M0c::cview_type M0cv;
            M0c m0c = m0;
            M0cv m0cv = m0c.cView();
            HouseholderMultEq_Helper<-2,V,M0cv,Mx,Vt>::call(
                u,beta,m0cv,mx,temp);
            NoAliasCopy(m0c,m0);
        }
    };

    // algo 90: Call inst
    template <class V, class M0, class Mx, class Vt>
    struct HouseholderMultEq_Helper<90,V,M0,Mx,Vt>
    {
        typedef typename V::real_type RT;
        static void call(const V& u, RT beta, M0& m0, Mx& mx, Vt& temp)
        { 
            InstHouseholderMultEq(
                u.xView(),beta,m0.xView(),mx.xView(),temp.xView()); 
        }
    };

    // algo 97: Conjugate mx
    template <class V, class M0, class Mx, class Vt>
    struct HouseholderMultEq_Helper<97,V,M0,Mx,Vt>
    {
        typedef typename V::real_type RT;
        static void call(const V& u, RT beta, M0& m0, Mx& mx, Vt& temp)
        { 
            typedef typename V::const_conjugate_type Vc;
            typedef typename M0::conjugate_type M0c;
            typedef typename Mx::conjugate_type Mxc;
            Vc uc = u.conjugate();
            M0c m0c = m0.conjugate();
            Mxc mxc = mx.conjugate();
            HouseholderMultEq_Helper<-2,Vc,M0c,Mxc,Vt>::call(
                uc,beta,m0c,mxc,temp);
        }
    };

    // algo -3: Select algorithm
    template <class V, class M0, class Mx, class Vt>
    struct HouseholderMultEq_Helper<-3,V,M0,Mx,Vt>
    {
        typedef typename V::real_type RT;
        static void call(const V& u, RT beta, M0& m0, Mx& mx, Vt& temp)
        { 
            const int algo = 11;
#ifdef PRINTALGO_HOUSE
            std::cout<<"Inline HouseholderMultEq\n";
            std::cout<<"u = "<<u<<std::endl;
            std::cout<<"beta = "<<beta<<std::endl;
            std::cout<<"m0 = "<<m0<<std::endl;
            std::cout<<"mx = "<<mx<<std::endl;
            std::cout<<"temp = "<<temp<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            HouseholderMultEq_Helper<algo,V,M0,Mx,Vt>::call(u,beta,m0,mx,temp);
#ifdef PRINTALGO_HOUSE
            std::cout<<"m0 => "<<m0<<std::endl;
            std::cout<<"mx => "<<mx<<std::endl;
#endif
        }
    };

    // algo -2: Check for inst
    template <class V, class M0, class Mx, class Vt>
    struct HouseholderMultEq_Helper<-2,V,M0,Mx,Vt>
    {
        typedef typename V::real_type RT;
        static void call(const V& u, RT beta, M0& m0, Mx& mx, Vt& temp)
        { 
            const bool inst = 
                (V::_size == TMV_UNKNOWN || V::_size > 16) &&
                (Mx::_colsize == TMV_UNKNOWN || Mx::_colsize > 16) &&
                (Mx::_rowsize == TMV_UNKNOWN || Mx::_rowsize > 16) &&
                (M0::_size == TMV_UNKNOWN || M0::_size > 16) &&
                Traits<RT>::isinst;
            const int algo = 
                Mx::_conj ? 97 : 
                inst ? (M0::_conj ? 82 : 90) :
                -3;
            HouseholderMultEq_Helper<algo,V,M0,Mx,Vt>::call(u,beta,m0,mx,temp);
        }
    };

    template <class V, class M0, class Mx, class Vt>
    struct HouseholderMultEq_Helper<-1,V,M0,Mx,Vt>
    {
        typedef typename V::real_type RT;
        static void call(const V& u, RT beta, M0& m0, Mx& mx, Vt& temp)
        { HouseholderMultEq_Helper<-2,V,M0,Mx,Vt>::call(u,beta,m0,mx,temp); }
    };

    template <class V, class M0, class Mx, class Vt>
    static inline void InlineHouseholderMultEq(
        const BaseVector_Calc<V>& u, typename V::real_type beta,
        BaseVector_Mutable<M0>& m0, BaseMatrix_Rec_Mutable<Mx>& mx,
        BaseVector_Mutable<Vt>& temp)
    {
        TMVStaticAssert((Sizes<V::_size,Mx::_colsize>::same));
        TMVStaticAssert((Sizes<M0::_size,Mx::_rowsize>::same));
        TMVStaticAssert((Sizes<Vt::_size,Mx::_rowsize>::same));
        TMVStaticAssert(Vt::_step == 1 || Vt::_step == TMV_UNKNOWN);
        TMVAssert(u.size() == mx.colsize());
        TMVAssert(m0.size() == mx.rowsize());
        TMVAssert(temp.size() == mx.rowsize());
        TMVAssert(temp.step() == 1);
        typedef typename V::const_cview_type Vv;
        typedef typename M0::cview_type M0v;
        typedef typename Mx::cview_type Mxv;
        typedef typename Vt::cview_type Vtv;
        TMV_MAYBE_CREF(V,Vv) uv = u.cView();
        TMV_MAYBE_REF(M0,M0v) m0v = m0.cView();
        TMV_MAYBE_REF(Mx,Mxv) mxv = mx.cView();
        TMV_MAYBE_REF(Vt,Vtv) tv = temp.cView();
        HouseholderMultEq_Helper<-3,Vv,M0v,Mxv,Vtv>::call(uv,beta,m0v,mxv,tv);
    }

    template <class V, class M0, class Mx, class Vt>
    static inline void HouseholderMultEq(
        const BaseVector_Calc<V>& u, typename V::real_type beta,
        BaseVector_Mutable<M0>& m0, BaseMatrix_Rec_Mutable<Mx>& mx,
        BaseVector_Mutable<Vt>& temp)
    {
        TMVStaticAssert((Sizes<V::_size,Mx::_colsize>::same));
        TMVStaticAssert((Sizes<M0::_size,Mx::_rowsize>::same));
        TMVStaticAssert((Sizes<Vt::_size,Mx::_rowsize>::same));
        TMVStaticAssert(Vt::_step == 1 || Vt::_step == TMV_UNKNOWN);
        TMVAssert(u.size() == mx.colsize());
        TMVAssert(m0.size() == mx.rowsize());
        TMVAssert(temp.size() == mx.rowsize());
        TMVAssert(temp.step() == 1);
        typedef typename V::const_cview_type Vv;
        typedef typename M0::cview_type M0v;
        typedef typename Mx::cview_type Mxv;
        typedef typename Vt::cview_type Vtv;
        TMV_MAYBE_CREF(V,Vv) uv = u.cView();
        TMV_MAYBE_REF(M0,M0v) m0v = m0.cView();
        TMV_MAYBE_REF(Mx,Mxv) mxv = mx.cView();
        TMV_MAYBE_REF(Vt,Vtv) tv = temp.cView();
        HouseholderMultEq_Helper<-2,Vv,M0v,Mxv,Vtv>::call(uv,beta,m0v,mxv,tv);
    }

    //
    // BlockHouseholderAugment
    // Augment a block householder (I - YZYt) matrix with one more 
    // householder matrix given by the last column of Y and beta.
    // This routine updates the last column of Z such that the latest
    // Householder matrix is multiplied on the right.
    //
    // Note: The convention I use for conjugation is that the non-conjugated
    // Z always corresponds to the non-conjugated Y.
    // So if Y is conj, you should use Z.conjugate() as its corresponding
    // Z matrix.
    //
    
    template <int algo, class M1, class M2>
    struct BlockHouseholderAugment_Helper;

    // algo 11: Normal case
    template <class M1, class M2>
    struct BlockHouseholderAugment_Helper<11,M1,M2>
    {
        typedef typename M2::real_type RT;
        static void call(const M1& Y, M2& Z, RT beta)
        {
            // All but the last columns of the input matrices, Y,Z are such 
            // that I - Y'Z'Y't is a Block Householder matrix (the product of 
            // several individual Householder matrices).
            // The last column of Y has the vector for the next Householder 
            // matrix.
            //
            // Z is output such that:
            // I - Y Z Yt = (I - Y' Z' Y't) (I - beta v vt)
            //
            // Y Z Yt = Y' Z' Y't + beta v vt - beta Y' Z' Y't v vt
            //
            // Blocking Y as [ Y' v ] and Z as [ Z'  z   ]
            //                                 [ 0  beta ]
            // Y Z Yt = [ Y' v ] [ Z' Y't + z vt ]
            //                   [    beta vt    ]
            //        = Y' Z' Y't + beta v vt + Y' z vt
            //
            // Comparing these equations, we find:
            // z = -beta Z' Y't v
            //
            // Remember that the first element of v is not stored, but rather
            // is assumed to be 1.

            int M = Y.colsize();
            int N = Y.rowsize()-1; // # of cols already computed
#ifdef XDEBUG_HOUSE
            typedef typename M1::value_type T;
            Matrix<T> Y0(Y);
            Y0.upperTri().setToIdentity();
            Matrix<T> Z0(Z);
            Matrix<T> H0 = T(1) - 
                Y0.colRange(0,N)*Z.subTriMatrix(0,N)*Y0.colRange(0,N).adjoint();
            Matrix<T> H1(Y.colsize(),Y.colsize());
            H1.setToIdentity();
            H1.subMatrix(N,M,N,M) -= beta * (Y0.col(N,N,M)^Y0.col(N,N,M).conjugate());
            Matrix<T> H2 = H0*H1;
#endif

            if (beta == RT(0)) {
                Z.col(N,0,N+1).setZero();
            } else {
                if (N > 0) {
                    typename M1::const_col_sub_type v = Y.col(N,N+1,M);
                    typename M2::col_sub_type z = Z.col(N,0,N);
                    z = Y.subMatrix(N+1,M,0,N).adjoint()*v;
                    z += Y.row(N,0,N).conjugate();
                    z = -beta * Z.subTriMatrix(0,N) * z;
                }
                Z.ref(N,N) = beta;
            }
#ifdef XDEBUG_HOUSE
            Matrix<T> YY(Y);
            YY.upperTri().setToIdentity();
            Matrix<T> HH = T(1) - YY * Z * YY.adjoint();
            if (!(Norm(HH-H2) <= 0.001*Norm(H2))) {
                cerr<<"BlockHouseholderAugment\n";
                cerr<<"Input: Y = "<<Y0<<endl;
                cerr<<"beta = "<<beta<<endl;
                cerr<<"Z = "<<Z0<<endl;
                cerr<<"H0 = "<<H0<<endl;
                cerr<<"H1 = "<<H1<<endl;
                cerr<<"H0*H1 = "<<H2<<endl;
                cerr<<"Output: Y = "<<Y<<endl;
                cerr<<"Z = "<<Z<<endl;
                cerr<<"H = "<<HH<<endl;
                abort();
            }
#endif
        }
    };

    // algo 90: Call inst
    template <class M1, class M2>
    struct BlockHouseholderAugment_Helper<90,M1,M2>
    {
        typedef typename M2::real_type RT;
        static void call(const M1& Y, M2& Z, RT beta)
        { InstBlockHouseholderAugment(Y.xView(),Z.xView(),beta); }
    };

    // algo 97: Conjugate Y,Z
    template <class M1, class M2>
    struct BlockHouseholderAugment_Helper<97,M1,M2>
    {
        typedef typename M2::real_type RT;
        static void call(const M1& Y, M2& Z, RT beta)
        { 
            typedef typename M1::const_nonconj_type M1c;
            typedef typename M2::nonconj_type M2c;
            M1c Yc = Y.nonConj();
            M2c Zc = Z.nonConj();
            BlockHouseholderAugment_Helper<-2,M1c,M2c>::call(Yc,Zc,beta);
        }
    };

    // algo -3: Select algorithm
    template <class M1, class M2>
    struct BlockHouseholderAugment_Helper<-3,M1,M2>
    {
        typedef typename M2::real_type RT;
        static void call(const M1& Y, M2& Z, RT beta)
        { 
            const int algo = 11;
#ifdef PRINTALGO_HOUSE
            std::cout<<"Inline BlockHouseholderAugment\n";
            std::cout<<"Y = "<<Y<<std::endl;
            std::cout<<"Z = "<<Z<<std::endl;
            std::cout<<"beta = "<<beta<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            BlockHouseholderAugment_Helper<algo,M1,M2>::call(Y,Z,beta);
#ifdef PRINTALGO_HOUSE
            std::cout<<"Z => "<<Z<<std::endl;
#endif
        }
    };

    // algo -2: Check for inst
    template <class M1, class M2>
    struct BlockHouseholderAugment_Helper<-2,M1,M2>
    {
        typedef typename M2::real_type RT;
        static void call(const M1& Y, M2& Z, RT beta)
        { 
            const bool inst = 
                (M1::_colsize == TMV_UNKNOWN || M1::_colsize > 16) &&
                (M1::_rowsize == TMV_UNKNOWN || M1::_rowsize > 16) &&
                (M2::_size == TMV_UNKNOWN || M2::_size > 16) &&
                Traits<RT>::isinst;
            const int algo = 
                M1::_conj || M2::_conj ? 97 : 
                inst ? 90 :
                -3;
            BlockHouseholderAugment_Helper<algo,M1,M2>::call(Y,Z,beta);
        }
    };

    template <class M1, class M2>
    struct BlockHouseholderAugment_Helper<-1,M1,M2>
    {
        typedef typename M2::real_type RT;
        static void call(const M1& Y, M2& Z, RT beta)
        { BlockHouseholderAugment_Helper<-2,M1,M2>::call(Y,Z,beta); }
    };

    template <class M1, class M2>
    static inline void InlineBlockHouseholderAugment(
        const BaseMatrix_Rec<M1>& Y, BaseMatrix_Tri_Mutable<M2>& Z,
        typename M2::real_type beta)
    {
        TMVStaticAssert((Sizes<M2::_size,M1::_rowsize>::same));
        TMVAssert(Z.size() == Y.rowsize());
        TMVAssert(Y.rowsize() > 0);
        TMVAssert(Y.colsize() > 0);
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) Yv = Y.cView();
        TMV_MAYBE_REF(M2,M2v) Zv = Z.cView();
        BlockHouseholderAugment_Helper<-3,M1v,M2v>::call(Yv,Zv,beta);
    }

    template <class M1, class M2>
    static inline void BlockHouseholderAugment(
        const BaseMatrix_Rec<M1>& Y, BaseMatrix_Tri_Mutable<M2>& Z,
        typename M2::real_type beta)
    {
        TMVStaticAssert((Sizes<M2::_size,M1::_rowsize>::same));
        TMVAssert(Z.size() == Y.rowsize());
        TMVAssert(Y.rowsize() > 0);
        TMVAssert(Y.colsize() > 0);
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) Yv = Y.cView();
        TMV_MAYBE_REF(M2,M2v) Zv = Z.cView();
        BlockHouseholderAugment_Helper<-2,M1v,M2v>::call(Yv,Zv,beta);
    }

    //
    // BlockHouseholderMakeZ
    // This routine calculates the Z component of the BlockHouseholder
    // formulation for Q.  Y contains the v's for the Householder matrices,
    // and beta contains the beta's.  
    // The output BlockHouseholder matrix I-YZYt is the product 
    // H0 H1 ... HN.
    //

    template <int algo, class M1, class M2, class V>
    struct BlockHouseholderMakeZ_Helper;

    // algo 21: Recursive calculation
    template <class M1, class M2, class V>
    struct BlockHouseholderMakeZ_Helper<21,M1,M2,V>
    {
        static void call(const M1& Y, M2& Z, const V& beta)
        {
            // This routine calculates the Z component of the BlockHouseholder
            // formulation for Q.  Y contains the v's for the Householder 
            // matrices, and beta contains the beta's.  
            //
            // The output BlockHouseholder matrix I-YZYt is the product 
            // H0 H1 ... HN.
            //
            // Note that the Y matrix is really just the unit lower trapezoidal 
            // component of the input Y.

            typedef typename M1::value_type T;
            typedef typename M1::real_type RT;

            const int M = Y.colsize();
            const int N = Y.rowsize();

#ifdef XDEBUG_HOUSE
            Matrix<T> Y0(Y);
            Y0.upperTri().setToIdentity();
            Matrix<T> Z0(Z);
            Vector<T> beta0 = beta;
            Matrix<T> Htot(M,M);
            Htot.setToIdentity();
            for(int i=0;i<N;i++) {
                Matrix<T> H(M,M);
                H.setToIdentity();
                H.subMatrix(i,M,i,M) -= beta(i) * 
                    (Y0.col(i,i,M) ^ Y0.col(i,i,M).conjugate());
                Htot *= H.adjoint();
            }
#endif

            if (N > 2) {
                int j1 = (N+1)/2;
                typename M1::const_colrange_type Y1 = Y.colRange(0,j1);
                typename M2::subtrimatrix_type Z1 = Z.subTriMatrix(0,j1);
                typename V::const_subvector_type beta1 = beta.subVector(0,j1);
                BlockHouseholderMakeZ(Y1,Z1,beta1);

                typename M1::const_submatrix_type Y2 = Y.subMatrix(j1,M,j1,N);
                typename M2::subtrimatrix_type Z2 = Z.subTriMatrix(j1,N);
                typename V::const_subvector_type beta2 = beta.subVector(j1,N);
                BlockHouseholderMakeZ(Y2,Z2,beta2);

                // (I-Y1 Z1 Y1t)(I-Y2 Z2 Y2t) = 
                // I - Y1 Z1 Y1t - Y2 Z2 Y2t + Y1 Z1 Y1t Y2 Z2 Y2t
                // Y = [ Y1 Y2 ]
                // Z = [ Z1 Z3 ]
                //     [ 0  Z2 ]
                // Y Z Yt = [ Y1 Y2 ] [ Z1 Y1t + Z3 Y2t ]
                //                    [      Z2 Y2t     ]
                //        = Y1 Z1 Y1t + Y2 Z2 Y2t + Y1 Z3 Y2t
                // So, Z3 = -Z1 Y1t Y2 Z2
                // Remember that the Y's are Unit Lower Trapezoidal, so do the 
                // rectangle and triangle parts separately.
                //
                typename M2::submatrix_type Z3 = Z.subMatrix(0,j1,j1,N);
                Z3 = Y1.rowRange(j1,N).adjoint() *
                    Y.subMatrix(j1,N,j1,N).unitLowerTri();
                Z3 += Y1.rowRange(N,M).adjoint() * Y.subMatrix(N,M,j1,N);
                Z3 = -Z1*Z3;
                Z3 *= Z2;
            } else if (N==2) {
                // Z1 = Z(0,0) = beta(0)
                // Z2 = Z(1,1) = beta(1)
                // Y1 = ( 1 Y(1,0) Y.col(0,2,M) )
                // Y2 = ( 0   1    Y.col(1,2,M) )
                // Z3 = Z(0,1) = -Z1 Y1t Y2 Z2
                //    = -b0 b1 Y1* Y2
                const RT b0 = beta.cref(0);
                const RT b1 = beta.cref(1);
                Z.ref(0,0) = b0;
                Z.ref(1,1) = b1;
                T temp = Y.col(0,2,M).conjugate()*Y.col(1,2,M);
                temp += TMV_CONJ(Y.cref(1,0));
                Z.ref(0,1) = -b0*b1*temp;
            } else { // N == 1
                // I - Y Z Yt = I - beta* v vt
                // Therefore:
                // Y.col(0) = v
                // Z(0,0) = beta
                Z.ref(0,0) = beta.cref(0);
            }
#ifdef XDEBUG_HOUSE
            Matrix<T> Y2(Y);
            Y2.upperTri().setToIdentity();
            Matrix<T> Hnet = T(1) - Y2*Z*Y2.adjoint();
            if (!(Norm(Htot-Hnet) <= 0.001*Norm(Htot))) {
                cerr<<"BlockHouseholderMakeZ\n";
                cerr<<"Input: Y = "<<Y0<<endl;
                cerr<<"beta = "<<beta<<endl;
                cerr<<"H = "<<Htot<<endl;
                cerr<<"Output: Y = "<<Y2<<endl;
                cerr<<"Z = "<<Z<<endl;
                cerr<<"H = "<<Hnet<<endl;
                cerr<<"Norm(H-Htot) = "<<Norm(Htot-Hnet)<<endl;
                abort();
            }
#endif
        }
    };

    // algo 90: Call inst
    template <class M1, class M2, class V>
    struct BlockHouseholderMakeZ_Helper<90,M1,M2,V>
    {
        static void call(const M1& Y, M2& Z, const V& beta)
        { InstBlockHouseholderMakeZ(Y.xView(),Z.xView(),beta.xView()); }
    };

    // algo 97: Conjugate Y,Z
    template <class M1, class M2, class V>
    struct BlockHouseholderMakeZ_Helper<97,M1,M2,V>
    {
        static void call(const M1& Y, M2& Z, const V& beta)
        { 
            typedef typename M1::const_nonconj_type M1c;
            typedef typename M2::nonconj_type M2c;
            M1c Yc = Y.nonConj();
            M2c Zc = Z.nonConj();
            BlockHouseholderMakeZ_Helper<-2,M1c,M2c,V>::call(Yc,Zc,beta);
        }
    };

    // algo -3: Select algorithm
    template <class M1, class M2, class V>
    struct BlockHouseholderMakeZ_Helper<-3,M1,M2,V>
    {
        static void call(const M1& Y, M2& Z, const V& beta)
        { 
            const int algo = 21;
#ifdef PRINTALGO_HOUSE
            std::cout<<"Inline BlockHouseholderMakeZ\n";
            std::cout<<"Y = "<<Y<<std::endl;
            std::cout<<"Z = "<<Z<<std::endl;
            std::cout<<"beta = "<<beta<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            BlockHouseholderMakeZ_Helper<algo,M1,M2,V>::call(Y,Z,beta);
#ifdef PRINTALGO_HOUSE
            std::cout<<"Z => "<<Z<<std::endl;
#endif
        }
    };

    // algo -2: Check for inst
    template <class M1, class M2, class V>
    struct BlockHouseholderMakeZ_Helper<-2,M1,M2,V>
    {
        static void call(const M1& Y, M2& Z, const V& beta)
        { 
            typedef typename V::real_type RT;
            const bool inst = 
                (V::_size == TMV_UNKNOWN || V::_size > 16) &&
                (M1::_colsize == TMV_UNKNOWN || M1::_colsize > 16) &&
                (M1::_rowsize == TMV_UNKNOWN || M1::_rowsize > 16) &&
                (M2::_size == TMV_UNKNOWN || M2::_size > 16) &&
                Traits<RT>::isinst;
            const int algo = 
                M1::_conj || M2::_conj ? 97 : 
                inst ? 90 :
                -3;
            BlockHouseholderMakeZ_Helper<algo,M1,M2,V>::call(Y,Z,beta);
        }
    };

    template <class M1, class M2, class V>
    struct BlockHouseholderMakeZ_Helper<-1,M1,M2,V>
    {
        static void call(const M1& Y, M2& Z, const V& beta)
        { BlockHouseholderMakeZ_Helper<-2,M1,M2,V>::call(Y,Z,beta); }
    };

    template <class M1, class M2, class V>
    static inline void InlineBlockHouseholderMakeZ(
        const BaseMatrix_Rec<M1>& Y, BaseMatrix_Tri_Mutable<M2>& Z, 
        const BaseVector_Calc<V>& beta)
    {
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_size>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,V::_size>::same));
        TMVAssert(Y.colsize() >= Y.rowsize());
        TMVAssert(Y.rowsize() > 0);
        TMVAssert(Y.colsize() > 0);
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(M1,M1v) Yv = Y.cView();
        TMV_MAYBE_REF(M2,M2v) Zv = Z.cView();
        TMV_MAYBE_CREF(V,Vv) bv = beta.cView();
        BlockHouseholderMakeZ_Helper<-3,M1v,M2v,Vv>::call(Yv,Zv,bv);
    }

    template <class M1, class M2, class V>
    static inline void BlockHouseholderMakeZ(
        const BaseMatrix_Rec<M1>& Y, BaseMatrix_Tri_Mutable<M2>& Z, 
        const BaseVector_Calc<V>& beta)
    {
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_size>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,V::_size>::same));
        TMVAssert(Y.colsize() >= Y.rowsize());
        TMVAssert(Y.rowsize() > 0);
        TMVAssert(Y.colsize() > 0);
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(M1,M1v) Yv = Y.cView();
        TMV_MAYBE_REF(M2,M2v) Zv = Z.cView();
        TMV_MAYBE_CREF(V,Vv) bv = beta.cView();
        BlockHouseholderMakeZ_Helper<-2,M1v,M2v,Vv>::call(Yv,Zv,bv);
    }


    //
    // BlockHouseholderLMult
    // The input Y,Z are such that (I - YZYt) is a Block Householder matrix.
    // The upper square portion of Y is taken to be unit lower triangular.
    // ie. the diagonal and upper triangular portion are not referenced.
    // The routine then finds m <- (I - YZYt) m
    //

    template <int algo, class M1, class M2, class M3, class M4>
    struct BlockHouseholderLMult_Helper;

    // algo 11: Normal case
    template <class M1, class M2, class M3, class M4>
    struct BlockHouseholderLMult_Helper<11,M1,M2,M3,M4>
    {
        static void call(const M1& Y, const M2& Z, M3& m, M4& temp)
        {
#ifdef XDEBUG_HOUSE
            typedef typename M1::value_type T;
            typedef typename M3::value_type T3;
            Matrix<T> Y0(Y);
            Y0.upperTri().setToIdentity();
            Matrix<T> Z0(Z);
            Matrix<T3> m0(m);
            Matrix<T> H = T(1) - Y0*Z0*Y0.adjoint();
            Matrix<T3> Hm = H*m0;
#endif

            int M = Y.colsize();
            int N = Y.rowsize();
            const int cs = Sizes<M1::_colsize,M3::_colsize>::size;
            const int rs = Sizes<M1::_rowsize,M2::_size>::size;
            const int xs = M3::_rowsize;
            const int cs1 = IntTraits2<cs,rs>::diff;
            const int Si1 = M1::_stepi;
            const int Sj1 = M1::_stepj;
            const int C1 = M1::_conj;
            const int Si3 = M3::_stepi;
            const int Sj3 = M3::_stepj;
            const int C3 = M3::_conj;
            typedef typename M1::value_type T1;
            typedef typename M3::value_type T3;
            typedef typename MViewHelper<T1,Rec,rs,rs,Si1,Sj1,C1>::ctype M1a;
            typedef typename MViewHelper<T1,Rec,cs1,rs,Si1,Sj1,C1>::ctype M1b;
            typedef typename MViewHelper<T3,Rec,rs,xs,Si3,Sj3,C3>::type M3a;
            typedef typename MViewHelper<T3,Rec,cs1,xs,Si3,Sj3,C3>::type M3b;

            M1a Ya = Y.rowRange(0,N);
            M1b Yb = Y.rowRange(N,M);
            M3a ma = m.rowRange(0,N);
            M3b mb = m.rowRange(N,M);

            // temp = ZYtm
            temp = Ya.unitLowerTri().adjoint() * ma;
            temp += Yb.adjoint() * mb;
            temp = -Z * temp;
            ma += Ya.unitLowerTri() * temp;
            mb += Yb * temp;

#ifdef XDEBUG_HOUSE
            if (!(Norm(Hm-m) <= 0.001*Norm(m0)*Norm(H))) {
                cerr<<"BlockHouseholderLMult\n";
                cerr<<"Input: Y = "<<Y0<<endl;
                cerr<<"Z = "<<Z0<<endl;
                cerr<<"m = "<<m0<<endl;
                cerr<<"H = "<<H<<endl;
                cerr<<"Ya = "<<Ya<<endl;
                cerr<<"Yb = "<<Yb<<endl;
                cerr<<"temp = "<<temp<<endl;
                cerr<<"Output: m = "<<m<<endl;
                cerr<<"Hm = "<<Hm<<endl;
                abort();
            }
#endif
        }
    };

    // algo 90: Call inst
    template <class M1, class M2, class M3, class M4>
    struct BlockHouseholderLMult_Helper<90,M1,M2,M3,M4>
    {
        static void call(const M1& Y, const M2& Z, M3& m, M4& temp)
        {
            InstBlockHouseholderLMult(
                Y.xView(),Z.xView(),m.xView(),temp.xView()); 
        }
    };

    // algo 95: Conjugate Z to match Y
    template <class M1, class M2, class M3, class M4>
    struct BlockHouseholderLMult_Helper<95,M1,M2,M3,M4>
    {
        static void call(const M1& Y, const M2& Z, M3& m, M4& temp)
        { 
            typedef typename M2::const_conjugate_type M2c;
            M2c Zc = Z.conjugate();
            BlockHouseholderLMult_Helper<-2,M1,M2c,M3,M4>::call(Y,Zc,m,temp);
        }
    };

    // algo 97: Conjugate m
    template <class M1, class M2, class M3, class M4>
    struct BlockHouseholderLMult_Helper<97,M1,M2,M3,M4>
    {
        static void call(const M1& Y, const M2& Z, M3& m, M4& temp)
        { 
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            typedef typename M3::conjugate_type M3c;
            M1c Yc = Y.conjugate();
            M2c Zc = Z.conjugate();
            M3c mc = m.conjugate();
            BlockHouseholderLMult_Helper<-2,M1c,M2c,M3c,M4>::call(
                Yc,Zc,mc,temp);
        }
    };

    // algo -3: Select algorithm
    template <class M1, class M2, class M3, class M4>
    struct BlockHouseholderLMult_Helper<-3,M1,M2,M3,M4>
    {
        static void call(const M1& Y, const M2& Z, M3& m, M4& temp)
        { 
            const int algo = 11;
#ifdef PRINTALGO_HOUSE
            std::cout<<"Inline BlockHouseholderLMult\n";
            std::cout<<"Y = "<<Y<<std::endl;
            std::cout<<"Z = "<<Z<<std::endl;
            std::cout<<"m = "<<m<<std::endl;
            std::cout<<"temp = "<<temp<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            BlockHouseholderLMult_Helper<algo,M1,M2,M3,M4>::call(Y,Z,m,temp);
#ifdef PRINTALGO_HOUSE
            std::cout<<"m => "<<m<<std::endl;
#endif
        }
    };

    // algo -2: Check for inst
    template <class M1, class M2, class M3, class M4>
    struct BlockHouseholderLMult_Helper<-2,M1,M2,M3,M4>
    {
        static void call(const M1& Y, const M2& Z, M3& m, M4& temp)
        { 
            typedef typename M3::real_type RT;
            const bool inst = 
                (M1::_colsize == TMV_UNKNOWN || M1::_colsize > 16) &&
                (M1::_rowsize == TMV_UNKNOWN || M1::_rowsize > 16) &&
                (M2::_size == TMV_UNKNOWN || M2::_size > 16) &&
                (M3::_colsize == TMV_UNKNOWN || M3::_colsize > 16) &&
                (M3::_rowsize == TMV_UNKNOWN || M3::_rowsize > 16) &&
                Traits<RT>::isinst;
            const int algo = 
                M3::_conj ? 97 : 
                M2::_conj != int(M1::_conj) ? 95 :
                inst ? 90 :
                -3;
            BlockHouseholderLMult_Helper<algo,M1,M2,M3,M4>::call(Y,Z,m,temp);
        }
    };

    template <class M1, class M2, class M3, class M4>
    struct BlockHouseholderLMult_Helper<-1,M1,M2,M3,M4>
    {
        static void call(const M1& Y, const M2& Z, M3& m, M4& temp)
        { BlockHouseholderLMult_Helper<-2,M1,M2,M3,M4>::call(Y,Z,m,temp); }
    };

    template <class M1, class M2, class M3, class M4>
    static inline void InlineBlockHouseholderLMult(
        const BaseMatrix_Rec<M1>& Y, const BaseMatrix_Tri<M2>& Z,
        BaseMatrix_Rec_Mutable<M3>& m, BaseMatrix_Rec_Mutable<M4>& temp)
    {
        TMVStaticAssert((Sizes<M2::_size,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M3::_colsize,M1::_colsize>::same));
        TMVStaticAssert((Sizes<M4::_colsize,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M4::_rowsize,M3::_rowsize>::same));
        TMVAssert(Z.size() == Y.rowsize());
        TMVAssert(Y.rowsize() > 0);
        TMVAssert(Y.colsize() > 0);
        TMVAssert(m.colsize() == Y.colsize());
        TMVAssert(temp.colsize() == Y.rowsize());
        TMVAssert(temp.rowsize() == m.rowsize());
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        typedef typename M4::cview_type M4v;
        TMV_MAYBE_CREF(M1,M1v) Yv = Y.cView();
        TMV_MAYBE_CREF(M2,M2v) Zv = Z.cView();
        TMV_MAYBE_REF(M3,M3v) mv = m.cView();
        TMV_MAYBE_REF(M4,M4v) tv = temp.cView();
        BlockHouseholderLMult_Helper<-3,M1v,M2v,M3v,M4v>::call(Yv,Zv,mv,tv);
    }

    template <class M1, class M2, class M3, class M4>
    static inline void BlockHouseholderLMult(
        const BaseMatrix_Rec<M1>& Y, const BaseMatrix_Tri<M2>& Z,
        BaseMatrix_Rec_Mutable<M3>& m, BaseMatrix_Rec_Mutable<M4>& temp)
    {
        TMVStaticAssert((Sizes<M2::_size,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M3::_colsize,M1::_colsize>::same));
        TMVStaticAssert((Sizes<M4::_colsize,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M4::_rowsize,M3::_rowsize>::same));
        TMVAssert(Z.size() == Y.rowsize());
        TMVAssert(Y.rowsize() > 0);
        TMVAssert(Y.colsize() > 0);
        TMVAssert(m.colsize() == Y.colsize());
        TMVAssert(temp.colsize() == Y.rowsize());
        TMVAssert(temp.rowsize() == m.rowsize());
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        typedef typename M4::cview_type M4v;
        TMV_MAYBE_CREF(M1,M1v) Yv = Y.cView();
        TMV_MAYBE_CREF(M2,M2v) Zv = Z.cView();
        TMV_MAYBE_REF(M3,M3v) mv = m.cView();
        TMV_MAYBE_REF(M4,M4v) tv = temp.cView();
        BlockHouseholderLMult_Helper<-2,M1v,M2v,M3v,M4v>::call(Yv,Zv,mv,tv);
    }

    //
    // BlockHouseholderLDiv
    // The input Y,Z are such that (I - YZYt) is a Block Householder matrix.
    // The upper square portion of Y is taken to be unit lower triangular.
    // ie. the diagonal and upper triangular portion are not referenced.
    // The routine then finds m <- (I - YZYt)^-1 m = (I - YZtYt) m
    //
    // TODO: Split m into ma and mb in calling routing like for Block2 version.
    // Then calling routing might be able to know one of the sizes and
    // use a SmallMatrixView instead of MatrixView.

    template <int algo, class M1, class M2, class M3, class M4>
    struct BlockHouseholderLDiv_Helper;

    // algo 11: Normal case
    template <class M1, class M2, class M3, class M4>
    struct BlockHouseholderLDiv_Helper<11,M1,M2,M3,M4>
    {
        static void call(const M1& Y, const M2& Z, M3& m, M4& temp)
        {
            // The routine finds m <- (I - YZYt)^-1 m
            // = (I - YZtYt) m

            TMVAssert(Z.size() == Y.rowsize());
            TMVAssert(Y.rowsize() > 0);
            TMVAssert(Y.colsize() > 0);
            TMVAssert(m.colsize() == Y.colsize());
            TMVAssert(temp.colsize() == Y.rowsize());
            TMVAssert(temp.rowsize() == m.rowsize());
#ifdef XDEBUG_HOUSE
            typedef typename M1::value_type T;
            typedef typename M3::value_type T3;
            Matrix<T> Y0(Y);
            Y0.upperTri().setToIdentity();
            Matrix<T> Z0(Z.adjoint());
            Matrix<T3> m0(m);
            Matrix<T> Hinv = T(1) - Y0*Z0*Y0.adjoint();
            Matrix<T3> Hm = Hinv*m0;
#endif

            int M = Y.colsize();
            int N = Y.rowsize();
            const int cs = Sizes<M1::_colsize,M3::_colsize>::size;
            const int rs = Sizes<M1::_rowsize,M2::_size>::size;
            const int xs = M3::_rowsize;
            const int cs1 = IntTraits2<cs,rs>::diff;
            const int Si1 = M1::_stepi;
            const int Sj1 = M1::_stepj;
            const int C1 = M1::_conj;
            const int Si3 = M3::_stepi;
            const int Sj3 = M3::_stepj;
            const int C3 = M3::_conj;
            typedef typename M1::value_type T1;
            typedef typename M3::value_type T3;
            typedef typename MViewHelper<T1,Rec,rs,rs,Si1,Sj1,C1>::ctype M1a;
            typedef typename MViewHelper<T1,Rec,cs1,rs,Si1,Sj1,C1>::ctype M1b;
            typedef typename MViewHelper<T3,Rec,rs,xs,Si3,Sj3,C3>::type M3a;
            typedef typename MViewHelper<T3,Rec,cs1,xs,Si3,Sj3,C3>::type M3b;

            M1a Ya = Y.rowRange(0,N);
            M1b Yb = Y.rowRange(N,M);
            M3a ma = m.rowRange(0,N);
            M3b mb = m.rowRange(N,M);

            // temp = ZtYtm
            temp = Ya.unitLowerTri().adjoint() * ma;
            temp += Yb.adjoint() * mb;
            temp = -Z.adjoint() * temp;
            ma += Ya.unitLowerTri() * temp;
            mb += Yb * temp;

#ifdef XDEBUG_HOUSE
            if (!(Norm(Hm-m) <= 0.001*Norm(m0)*Norm(Hinv))) {
                cerr<<"BlockHouseholderLDiv\n";
                cerr<<"Input: Y = "<<Y0<<endl;
                cerr<<"Z = "<<Z0<<endl;
                cerr<<"m = "<<m0<<endl;
                cerr<<"Hinv = "<<Hinv<<endl;
                cerr<<"Output: m = "<<m<<endl;
                abort();
            }
#endif
        }
    };

    // algo 90: Call inst
    template <class M1, class M2, class M3, class M4>
    struct BlockHouseholderLDiv_Helper<90,M1,M2,M3,M4>
    {
        static void call(const M1& Y, const M2& Z, M3& m, M4& temp)
        {
            InstBlockHouseholderLDiv(
                Y.xView(),Z.xView(),m.xView(),temp.xView()); 
        }
    };

    // algo 95: Conjugate Z to match Y
    template <class M1, class M2, class M3, class M4>
    struct BlockHouseholderLDiv_Helper<95,M1,M2,M3,M4>
    {
        static void call(const M1& Y, const M2& Z, M3& m, M4& temp)
        { 
            typedef typename M2::const_conjugate_type M2c;
            M2c Zc = Z.conjugate();
            BlockHouseholderLDiv_Helper<-2,M1,M2c,M3,M4>::call(Y,Zc,m,temp);
        }
    };

    // algo 97: Conjugate m
    template <class M1, class M2, class M3, class M4>
    struct BlockHouseholderLDiv_Helper<97,M1,M2,M3,M4>
    {
        static void call(const M1& Y, const M2& Z, M3& m, M4& temp)
        { 
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            typedef typename M3::conjugate_type M3c;
            M1c Yc = Y.conjugate();
            M2c Zc = Z.conjugate();
            M3c mc = m.conjugate();
            BlockHouseholderLDiv_Helper<-2,M1c,M2c,M3c,M4>::call(
                Yc,Zc,mc,temp);
        }
    };

    // algo -3: Select algorithm
    template <class M1, class M2, class M3, class M4>
    struct BlockHouseholderLDiv_Helper<-3,M1,M2,M3,M4>
    {
        static void call(const M1& Y, const M2& Z, M3& m, M4& temp)
        { 
            const int algo = 11;
#ifdef PRINTALGO_HOUSE
            std::cout<<"Inline BlockHouseholderLDiv\n";
            std::cout<<"Y = "<<Y<<std::endl;
            std::cout<<"Z = "<<Z<<std::endl;
            std::cout<<"m = "<<m<<std::endl;
            std::cout<<"temp = "<<temp<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            BlockHouseholderLDiv_Helper<algo,M1,M2,M3,M4>::call(Y,Z,m,temp);
#ifdef PRINTALGO_HOUSE
            std::cout<<"m => "<<m<<std::endl;
#endif
        }
    };

    // algo -2: Check for inst
    template <class M1, class M2, class M3, class M4>
    struct BlockHouseholderLDiv_Helper<-2,M1,M2,M3,M4>
    {
        static void call(const M1& Y, const M2& Z, M3& m, M4& temp)
        { 
            typedef typename M3::real_type RT;
            const bool inst = 
                (M1::_colsize == TMV_UNKNOWN || M1::_colsize > 16) &&
                (M1::_rowsize == TMV_UNKNOWN || M1::_rowsize > 16) &&
                (M2::_size == TMV_UNKNOWN || M2::_size > 16) &&
                (M3::_colsize == TMV_UNKNOWN || M3::_colsize > 16) &&
                (M3::_rowsize == TMV_UNKNOWN || M3::_rowsize > 16) &&
                Traits<RT>::isinst;
            const int algo = 
                M3::_conj ? 97 : 
                M2::_conj != int(M1::_conj) ? 95 :
                inst ? 90 :
                -3;
            BlockHouseholderLDiv_Helper<algo,M1,M2,M3,M4>::call(Y,Z,m,temp);
        }
    };

    template <class M1, class M2, class M3, class M4>
    struct BlockHouseholderLDiv_Helper<-1,M1,M2,M3,M4>
    {
        static void call(const M1& Y, const M2& Z, M3& m, M4& temp)
        { BlockHouseholderLDiv_Helper<-2,M1,M2,M3,M4>::call(Y,Z,m,temp); }
    };

    template <class M1, class M2, class M3, class M4>
    static inline void InlineBlockHouseholderLDiv(
        const BaseMatrix_Rec<M1>& Y, const BaseMatrix_Tri<M2>& Z,
        BaseMatrix_Rec_Mutable<M3>& m, BaseMatrix_Rec_Mutable<M4>& temp)
    {
        TMVStaticAssert((Sizes<M2::_size,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M3::_colsize,M1::_colsize>::same));
        TMVStaticAssert((Sizes<M4::_colsize,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M4::_rowsize,M3::_rowsize>::same));
        TMVAssert(Z.size() == Y.rowsize());
        TMVAssert(Y.rowsize() > 0);
        TMVAssert(Y.colsize() > 0);
        TMVAssert(m.colsize() == Y.colsize());
        TMVAssert(temp.colsize() == Y.rowsize());
        TMVAssert(temp.rowsize() == m.rowsize());
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        typedef typename M4::cview_type M4v;
        TMV_MAYBE_CREF(M1,M1v) Yv = Y.cView();
        TMV_MAYBE_CREF(M2,M2v) Zv = Z.cView();
        TMV_MAYBE_REF(M3,M3v) mv = m.cView();
        TMV_MAYBE_REF(M4,M4v) tv = temp.cView();
        BlockHouseholderLDiv_Helper<-3,M1v,M2v,M3v,M4v>::call(Yv,Zv,mv,tv);
    }

    template <class M1, class M2, class M3, class M4>
    static inline void BlockHouseholderLDiv(
        const BaseMatrix_Rec<M1>& Y, const BaseMatrix_Tri<M2>& Z,
        BaseMatrix_Rec_Mutable<M3>& m, BaseMatrix_Rec_Mutable<M4>& temp)
    {
        TMVStaticAssert((Sizes<M2::_size,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M3::_colsize,M1::_colsize>::same));
        TMVStaticAssert((Sizes<M4::_colsize,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M4::_rowsize,M3::_rowsize>::same));
        TMVAssert(Z.size() == Y.rowsize());
        TMVAssert(Y.rowsize() > 0);
        TMVAssert(Y.colsize() > 0);
        TMVAssert(m.colsize() == Y.colsize());
        TMVAssert(temp.colsize() == Y.rowsize());
        TMVAssert(temp.rowsize() == m.rowsize());
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        typedef typename M4::cview_type M4v;
        TMV_MAYBE_CREF(M1,M1v) Yv = Y.cView();
        TMV_MAYBE_CREF(M2,M2v) Zv = Z.cView();
        TMV_MAYBE_REF(M3,M3v) mv = m.cView();
        TMV_MAYBE_REF(M4,M4v) tv = temp.cView();
        BlockHouseholderLDiv_Helper<-2,M1v,M2v,M3v,M4v>::call(Yv,Zv,mv,tv);
    }


    //
    // BlockHouseholderUnpack
    // The BlockHouseholder matrix H is defined by Y,Z.
    // Y is unpacked in place.
    //

    template <int algo, class M1, class M2, class M3>
    struct BlockHouseholderUnpack_Helper;

    // algo 11: Normal case
    template <class M1, class M2, class M3>
    struct BlockHouseholderUnpack_Helper<11,M1,M2,M3>
    {
        static void call(M1& Y, const M2& Z, M3& temp)
        {
            // The BlockHouseholder matrix H is defined by Y,Z.
            // Y is unpacked in place.
            TMVAssert(Y.colsize() > 0);
            TMVAssert(Y.rowsize() > 0);
            TMVAssert(Y.colsize() >= Y.rowsize());
            TMVAssert(Y.rowsize() == Z.size());
            TMVStaticAssert(M3::_upper);
            TMVStaticAssert(!M3::_unit);
            TMVStaticAssert(!M3::_unknowndiag);
            TMVAssert(temp.size() == Y.rowsize());

            int M = Y.colsize();
            int N = Y.rowsize();
            const int cs = M1::_colsize;
            const int rs = Sizes<M1::_rowsize,M2::_size>::size;
            const int cs1 = IntTraits2<cs,rs>::diff;
            const int Si1 = M1::_stepi;
            const int Sj1 = M1::_stepj;
            const int C1 = M1::_conj;
            typedef typename M1::value_type T;
            typedef typename MViewHelper<T,Rec,rs,rs,Si1,Sj1,C1>::type M1a;
            typedef typename MViewHelper<T,Rec,cs1,rs,Si1,Sj1,C1>::type M1b;

            M1a Ya = Y.rowRange(0,N);
            M1b Yb = Y.rowRange(N,M);

            // Make Y equal to:
            // Ht [ I ] = (I - YZYt) [ I ]
            //    [ 0 ]              [ 0 ]
            temp = -Z * Ya.unitLowerTri().adjoint();
            Yb *= temp;
            Ya = Ya.unitLowerTri() * temp;
            Ya.diag().addToAll(T(1));
        }
    };

    // algo 90: Call inst
    template <class M1, class M2, class M3>
    struct BlockHouseholderUnpack_Helper<90,M1,M2,M3>
    {
        static void call(M1& Y, const M2& Z, M3& temp)
        { InstBlockHouseholderUnpack(Y.xView(),Z.xView(),temp.xView()); }
    };

    // algo 97: Conjugate Y,Z
    template <class M1, class M2, class M3>
    struct BlockHouseholderUnpack_Helper<97,M1,M2,M3>
    {
        static void call(M1& Y, const M2& Z, M3& temp)
        { 
            typedef typename M1::nonconj_type M1c;
            typedef typename M2::const_nonconj_type M2c;
            M1c Yc = Y.nonConj();
            M2c Zc = Z.nonConj();
            BlockHouseholderUnpack_Helper<-2,M1c,M2c,M3>::call(Yc,Zc,temp);
        }
    };

    // algo -3: Select algorithm
    template <class M1, class M2, class M3>
    struct BlockHouseholderUnpack_Helper<-3,M1,M2,M3>
    {
        static void call(M1& Y, const M2& Z, M3& temp)
        { 
            const int algo = 11;
#ifdef PRINTALGO_HOUSE
            std::cout<<"Inline BlockHouseholderUnPack\n";
            std::cout<<"Y = "<<Y<<std::endl;
            std::cout<<"Z = "<<Z<<std::endl;
            std::cout<<"temp = "<<temp<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            BlockHouseholderUnpack_Helper<algo,M1,M2,M3>::call(Y,Z,temp);
#ifdef PRINTALGO_HOUSE
            std::cout<<"Y => "<<Y<<std::endl;
#endif
        }
    };

    // algo -2: Check for inst
    template <class M1, class M2, class M3>
    struct BlockHouseholderUnpack_Helper<-2,M1,M2,M3>
    {
        static void call(M1& Y, const M2& Z, M3& temp)
        { 
            typedef typename M2::real_type RT;
            const bool inst = 
                (M1::_colsize == TMV_UNKNOWN || M1::_colsize > 16) &&
                (M1::_rowsize == TMV_UNKNOWN || M1::_rowsize > 16) &&
                (M2::_size == TMV_UNKNOWN || M2::_size > 16) &&
                (M3::_size == TMV_UNKNOWN || M3::_size > 16) &&
                Traits<RT>::isinst;
            const int algo = 
                M1::_conj || M2::_conj ? 97 : 
                inst ? 90 :
                -3;
            BlockHouseholderUnpack_Helper<algo,M1,M2,M3>::call(Y,Z,temp);
        }
    };

    template <class M1, class M2, class M3>
    struct BlockHouseholderUnpack_Helper<-1,M1,M2,M3>
    {
        static void call(M1& Y, const M2& Z, M3& temp)
        { BlockHouseholderUnpack_Helper<-2,M1,M2,M3>::call(Y,Z,temp); }
    };

    template <class M1, class M2, class M3>
    static inline void InlineBlockHouseholderUnpack(
        BaseMatrix_Rec_Mutable<M1>& Y, const BaseMatrix_Tri<M2>& Z,
        BaseMatrix_Tri_Mutable<M3>& temp)
    {
        TMVStaticAssert(M3::_upper);
        TMVStaticAssert(!M3::_unit);
        TMVStaticAssert(!M3::_unknowndiag);
        TMVStaticAssert((Sizes<M2::_size,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M3::_size,M1::_rowsize>::same));
        TMVAssert(Z.size() == Y.rowsize());
        TMVAssert(Y.rowsize() > 0);
        TMVAssert(Y.colsize() > 0);
        TMVAssert(Y.colsize() >= Y.rowsize());
        TMVAssert(temp.size() == Y.rowsize());
        typedef typename M1::cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_REF(M1,M1v) Yv = Y.cView();
        TMV_MAYBE_CREF(M2,M2v) Zv = Z.cView();
        TMV_MAYBE_REF(M3,M3v) tv = temp.cView();
        BlockHouseholderUnpack_Helper<-3,M1v,M2v,M3v>::call(Yv,Zv,tv);
    }

    template <class M1, class M2, class M3>
    static inline void BlockHouseholderUnpack(
        BaseMatrix_Rec_Mutable<M1>& Y, const BaseMatrix_Tri<M2>& Z,
        BaseMatrix_Tri_Mutable<M3>& temp)
    {
        TMVStaticAssert(M3::_upper);
        TMVStaticAssert(!M3::_unit);
        TMVStaticAssert(!M3::_unknowndiag);
        TMVStaticAssert((Sizes<M2::_size,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M3::_size,M1::_rowsize>::same));
        TMVAssert(Z.size() == Y.rowsize());
        TMVAssert(Y.rowsize() > 0);
        TMVAssert(Y.colsize() > 0);
        TMVAssert(Y.colsize() >= Y.rowsize());
        TMVAssert(temp.size() == Y.rowsize());
        typedef typename M1::cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_REF(M1,M1v) Yv = Y.cView();
        TMV_MAYBE_CREF(M2,M2v) Zv = Z.cView();
        TMV_MAYBE_REF(M3,M3v) tv = temp.cView();
        BlockHouseholderUnpack_Helper<-2,M1v,M2v,M3v>::call(Yv,Zv,tv);
    }


    // The next set of functions that use Block2Householder...
    // are for QRUpdate and QRDowndate.
    // All of the Householder vectors are basically (1,0,0,0,0...,v)
    // where there are some arbitrary number of 0's between the initial
    // 1 and the real information in v.
    // Furthermore, when these are combined into a Y matrix, the 
    // v's are all the same size, and the upper part is just an identity
    // matrix.
    // This contrasts with the normal BlockHouseholder formulation with
    // the upper part of Y being unit-lower-triangular.
    
    //
    // Block2HouseholderAugment
    //
    
    template <int algo, class M1, class M2>
    struct Block2HouseholderAugment_Helper;

    // algo 11: Normal case
    template <class M1, class M2>
    struct Block2HouseholderAugment_Helper<11,M1,M2>
    {
        typedef typename M2::real_type RT;
        static void call(const M1& Y, M2& Z, RT beta)
        {
            // z = -beta Z' Y't v

            int N = Y.rowsize()-1; // # of cols already computed
#ifdef XDEBUG_HOUSE
            int M = Y.colsize();
            typedef typename M1::value_type T;
            Matrix<T> Y0(M+N+1,N+1);
            Y0.rowRange(0,N+1).setToIdentity();
            Y0.rowRange(N+1,M+N+1) = Y;
            Matrix<T> Z0(Z);
            Matrix<T> H0 = T(1) - 
                Y0.colRange(0,N)*Z.subTriMatrix(0,N)*Y0.colRange(0,N).adjoint();
            Matrix<T> H1(M,M);
            H1.setToIdentity();
            H1.subMatrix(N,M,N,M) -= beta * (Y0.col(N,N,M)^Y0.col(N,N,M).conjugate());
            Matrix<T> H2 = H0*H1;
#endif

            if (beta == RT(0)) {
                Z.col(N,0,N+1).setZero();
            } else {
                if (N > 0) {
                    typename M1::const_col_type v = Y.col(N);
                    typename M2::col_sub_type z = Z.col(N,0,N);
                    z = Y.colRange(0,N).adjoint()*v;
                    z = -beta * Z.subTriMatrix(0,N) * z;
                }
                Z.ref(N,N) = beta;
            }
#ifdef XDEBUG_HOUSE
            Matrix<T> YY(M+N+1,N+1);
            YY.rowRange(0,N+1).setToIdentity();
            YY.rowRange(N+1,M+N+1) = Y;
            Matrix<T> HH = T(1) - YY * Z * YY.adjoint();
            if (!(Norm(HH-H2) <= 0.001*Norm(H2))) {
                cerr<<"Block2HouseholderAugment\n";
                cerr<<"Input: Y = "<<Y0<<endl;
                cerr<<"beta = "<<beta<<endl;
                cerr<<"Z = "<<Z0<<endl;
                cerr<<"H0 = "<<H0<<endl;
                cerr<<"H1 = "<<H1<<endl;
                cerr<<"H0*H1 = "<<H2<<endl;
                cerr<<"Output: Y = "<<Y<<endl;
                cerr<<"Z = "<<Z<<endl;
                cerr<<"H = "<<HH<<endl;
                abort();
            }
#endif
        }
    };

    // algo 90: Call inst
    template <class M1, class M2>
    struct Block2HouseholderAugment_Helper<90,M1,M2>
    {
        typedef typename M2::real_type RT;
        static void call(const M1& Y, M2& Z, RT beta)
        { InstBlock2HouseholderAugment(Y.xView(),Z.xView(),beta); }
    };

    // algo 97: Conjugate Y,Z
    template <class M1, class M2>
    struct Block2HouseholderAugment_Helper<97,M1,M2>
    {
        typedef typename M2::real_type RT;
        static void call(const M1& Y, M2& Z, RT beta)
        { 
            typedef typename M1::const_nonconj_type M1c;
            typedef typename M2::nonconj_type M2c;
            M1c Yc = Y.nonConj();
            M2c Zc = Z.nonConj();
            Block2HouseholderAugment_Helper<-2,M1c,M2c>::call(Yc,Zc,beta);
        }
    };

    // algo -3: Select algorithm
    template <class M1, class M2>
    struct Block2HouseholderAugment_Helper<-3,M1,M2>
    {
        typedef typename M2::real_type RT;
        static void call(const M1& Y, M2& Z, RT beta)
        { 
            const int algo = 11;
#ifdef PRINTALGO_HOUSE
            std::cout<<"Inline Block2HouseholderAugment\n";
            std::cout<<"Y = "<<Y<<std::endl;
            std::cout<<"Z = "<<Z<<std::endl;
            std::cout<<"beta = "<<beta<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            Block2HouseholderAugment_Helper<algo,M1,M2>::call(Y,Z,beta);
#ifdef PRINTALGO_HOUSE
            std::cout<<"Z => "<<Z<<std::endl;
#endif
        }
    };

    // algo -2: Check for inst
    template <class M1, class M2>
    struct Block2HouseholderAugment_Helper<-2,M1,M2>
    {
        typedef typename M2::real_type RT;
        static void call(const M1& Y, M2& Z, RT beta)
        { 
            const bool inst = 
                (M1::_colsize == TMV_UNKNOWN || M1::_colsize > 16) &&
                (M1::_rowsize == TMV_UNKNOWN || M1::_rowsize > 16) &&
                (M2::_size == TMV_UNKNOWN || M2::_size > 16) &&
                Traits<RT>::isinst;
            const int algo = 
                M1::_conj || M2::_conj ? 97 : 
                inst ? 90 :
                -3;
            Block2HouseholderAugment_Helper<algo,M1,M2>::call(Y,Z,beta);
        }
    };

    template <class M1, class M2>
    struct Block2HouseholderAugment_Helper<-1,M1,M2>
    {
        typedef typename M2::real_type RT;
        static void call(const M1& Y, M2& Z, RT beta)
        { Block2HouseholderAugment_Helper<-2,M1,M2>::call(Y,Z,beta); }
    };

    template <class M1, class M2>
    static inline void InlineBlock2HouseholderAugment(
        const BaseMatrix_Rec<M1>& Y, BaseMatrix_Tri_Mutable<M2>& Z,
        typename M2::real_type beta)
    {
        TMVStaticAssert((Sizes<M2::_size,M1::_rowsize>::same));
        TMVAssert(Z.size() == Y.rowsize());
        TMVAssert(Y.rowsize() > 0);
        TMVAssert(Y.colsize() > 0);
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) Yv = Y.cView();
        TMV_MAYBE_REF(M2,M2v) Zv = Z.cView();
        Block2HouseholderAugment_Helper<-3,M1v,M2v>::call(Yv,Zv,beta);
    }

    template <class M1, class M2>
    static inline void Block2HouseholderAugment(
        const BaseMatrix_Rec<M1>& Y, BaseMatrix_Tri_Mutable<M2>& Z,
        typename M2::real_type beta)
    {
        TMVStaticAssert((Sizes<M2::_size,M1::_rowsize>::same));
        TMVAssert(Z.size() == Y.rowsize());
        TMVAssert(Y.rowsize() > 0);
        TMVAssert(Y.colsize() > 0);
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) Yv = Y.cView();
        TMV_MAYBE_REF(M2,M2v) Zv = Z.cView();
        Block2HouseholderAugment_Helper<-2,M1v,M2v>::call(Yv,Zv,beta);
    }

    //
    // Block2HouseholderMakeZ
    //

    template <int algo, class M1, class M2, class V>
    struct Block2HouseholderMakeZ_Helper;

    // algo 21: Recursive calculation
    template <class M1, class M2, class V>
    struct Block2HouseholderMakeZ_Helper<21,M1,M2,V>
    {
        static void call(const M1& Y, M2& Z, const V& beta)
        {
            typedef typename M1::value_type T;
            typedef typename M1::real_type RT;
            const int N = Y.rowsize();
#ifdef XDEBUG_HOUSE
            const int M = Y.colsize();
            Matrix<T> Y0(M+N,N);
            Y0.rowRange(0,N).setToIdentity();
            Y0.rowRange(N,M+N) = Y;
            Matrix<T> Z0(Z);
            Vector<T> beta0 = beta;
            Matrix<T> Htot(M+N,M+N);
            Htot.setToIdentity();
            for(int i=0;i<N;i++) {
                Matrix<T> H(M+N,M+N);
                H.setToIdentity();
                H.subMatrix(i,M+N,i,M+N) -= beta(i) * 
                    (Y0.col(i,i,M+N) ^ Y0.col(i,i,M+N).conjugate());
                Htot *= H.adjoint();
            }
#endif

            if (N > 2) {
                int j1 = (N+1)/2;
                typename M1::const_colrange_type Y1 = Y.colRange(0,j1);
                typename M2::subtrimatrix_type Z1 = Z.subTriMatrix(0,j1);
                typename V::const_subvector_type beta1 = beta.subVector(0,j1);
                Block2HouseholderMakeZ(Y1,Z1,beta1);

                typename M1::const_submatrix_type Y2 = Y.colRange(j1,N);
                typename M2::subtrimatrix_type Z2 = Z.subTriMatrix(j1,N);
                typename V::const_subvector_type beta2 = beta.subVector(j1,N);
                Block2HouseholderMakeZ(Y2,Z2,beta2);

                // Z3 = -Z1 Y1t Y2 Z2
                typename M2::submatrix_type Z3 = Z.subMatrix(0,j1,j1,N);
                Z3 = Y1.adjoint() * Y2;
                Z3 = -Z1*Z3;
                Z3 *= Z2;
            } else if (N==2) {
                // Z3 = -b0 b1 Y1* Y2
                const RT b0 = beta.cref(0);
                const RT b1 = beta.cref(1);
                Z.ref(0,0) = b0;
                Z.ref(1,1) = b1;
                T temp = Y.col(0).conjugate()*Y.col(1);
                Z.ref(0,1) = -b0*b1*temp;
            } else { // N == 1
                Z.ref(0,0) = beta.cref(0);
            }
#ifdef XDEBUG_HOUSE
            Matrix<T> Y2(M+N,N);
            Y2.rowRange(0,N).setToIdentity();
            Y2.rowRange(N,M+N) = Y;
            Matrix<T> Hnet = T(1) - Y2*Z*Y2.adjoint();
            if (!(Norm(Htot-Hnet) <= 0.001*Norm(Htot))) {
                cerr<<"Block2HouseholderMakeZ\n";
                cerr<<"Input: Y = "<<Y0<<endl;
                cerr<<"beta = "<<beta<<endl;
                cerr<<"H = "<<Htot<<endl;
                cerr<<"Output: Y = "<<Y2<<endl;
                cerr<<"Z = "<<Z<<endl;
                cerr<<"H = "<<Hnet<<endl;
                cerr<<"Norm(H-Htot) = "<<Norm(Htot-Hnet)<<endl;
                abort();
            }
#endif
        }
    };

    // algo 90: Call inst
    template <class M1, class M2, class V>
    struct Block2HouseholderMakeZ_Helper<90,M1,M2,V>
    {
        static void call(const M1& Y, M2& Z, const V& beta)
        { InstBlock2HouseholderMakeZ(Y.xView(),Z.xView(),beta.xView()); }
    };

    // algo 97: Conjugate Y,Z
    template <class M1, class M2, class V>
    struct Block2HouseholderMakeZ_Helper<97,M1,M2,V>
    {
        static void call(const M1& Y, M2& Z, const V& beta)
        { 
            typedef typename M1::const_nonconj_type M1c;
            typedef typename M2::nonconj_type M2c;
            M1c Yc = Y.nonConj();
            M2c Zc = Z.nonConj();
            Block2HouseholderMakeZ_Helper<-2,M1c,M2c,V>::call(Yc,Zc,beta);
        }
    };

    // algo -3: Select algorithm
    template <class M1, class M2, class V>
    struct Block2HouseholderMakeZ_Helper<-3,M1,M2,V>
    {
        static void call(const M1& Y, M2& Z, const V& beta)
        { 
            const int algo = 21;
#ifdef PRINTALGO_HOUSE
            std::cout<<"Inline Block2HouseholderMakeZ\n";
            std::cout<<"Y = "<<Y<<std::endl;
            std::cout<<"Z = "<<Z<<std::endl;
            std::cout<<"beta = "<<beta<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            Block2HouseholderMakeZ_Helper<algo,M1,M2,V>::call(Y,Z,beta);
#ifdef PRINTALGO_HOUSE
            std::cout<<"Z => "<<Z<<std::endl;
#endif
        }
    };

    // algo -2: Check for inst
    template <class M1, class M2, class V>
    struct Block2HouseholderMakeZ_Helper<-2,M1,M2,V>
    {
        static void call(const M1& Y, M2& Z, const V& beta)
        { 
            typedef typename V::real_type RT;
            const bool inst = 
                (V::_size == TMV_UNKNOWN || V::_size > 16) &&
                (M1::_colsize == TMV_UNKNOWN || M1::_colsize > 16) &&
                (M1::_rowsize == TMV_UNKNOWN || M1::_rowsize > 16) &&
                (M2::_size == TMV_UNKNOWN || M2::_size > 16) &&
                Traits<RT>::isinst;
            const int algo = 
                M1::_conj || M2::_conj ? 97 : 
                inst ? 90 :
                -3;
            Block2HouseholderMakeZ_Helper<algo,M1,M2,V>::call(Y,Z,beta);
        }
    };

    template <class M1, class M2, class V>
    struct Block2HouseholderMakeZ_Helper<-1,M1,M2,V>
    {
        static void call(const M1& Y, M2& Z, const V& beta)
        { Block2HouseholderMakeZ_Helper<-2,M1,M2,V>::call(Y,Z,beta); }
    };

    template <class M1, class M2, class V>
    static inline void InlineBlock2HouseholderMakeZ(
        const BaseMatrix_Rec<M1>& Y, BaseMatrix_Tri_Mutable<M2>& Z, 
        const BaseVector_Calc<V>& beta)
    {
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_size>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,V::_size>::same));
        TMVAssert(Y.rowsize() > 0);
        TMVAssert(Y.colsize() > 0);
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(M1,M1v) Yv = Y.cView();
        TMV_MAYBE_REF(M2,M2v) Zv = Z.cView();
        TMV_MAYBE_CREF(V,Vv) bv = beta.cView();
        Block2HouseholderMakeZ_Helper<-3,M1v,M2v,Vv>::call(Yv,Zv,bv);
    }

    template <class M1, class M2, class V>
    static inline void Block2HouseholderMakeZ(
        const BaseMatrix_Rec<M1>& Y, BaseMatrix_Tri_Mutable<M2>& Z, 
        const BaseVector_Calc<V>& beta)
    {
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_size>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,V::_size>::same));
        TMVAssert(Y.rowsize() > 0);
        TMVAssert(Y.colsize() > 0);
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(M1,M1v) Yv = Y.cView();
        TMV_MAYBE_REF(M2,M2v) Zv = Z.cView();
        TMV_MAYBE_CREF(V,Vv) bv = beta.cView();
        Block2HouseholderMakeZ_Helper<-2,M1v,M2v,Vv>::call(Yv,Zv,bv);
    }


    //
    // Block2HouseholderLMult
    //

    template <int algo, class M1, class M2, class M3, class M4, class M5>
    struct Block2HouseholderLMult_Helper;

    // algo 11: Normal case
    template <class M1, class M2, class M3, class M4, class M5>
    struct Block2HouseholderLMult_Helper<11,M1,M2,M3,M4,M5>
    {
        static void call(const M1& Y, const M2& Z, M3& ma, M4& mb, M5& temp)
        {
#ifdef XDEBUG_HOUSE
            int M = Y.colsize();
            int N = Y.rowsize();
            int K = ma.rowsize();

            typedef typename M1::value_type T;
            typedef typename M3::value_type T3;
            Matrix<T> Y0(M+N,N);
            Y0.rowRange(0,N).setToIdentity();
            Y0.rowRange(N,M+N) = Y;
            Matrix<T> Z0(Z);
            Matrix<T3> m0(M+N,K);
            m0.rowRange(0,N) = ma;
            m0.rowRange(N,N+M) = mb;
            Matrix<T> H = T(1) - Y0*Z0*Y0.adjoint();
            Matrix<T3> Hm = H*m0;
#endif

            // temp = -ZYtm
            // Remember that the full Y includes an identity matrix on top.
            temp = Y.adjoint() * mb;
            temp += ma;
            temp = -Z * temp;
            ma += temp;
            mb += Y * temp;

#ifdef XDEBUG_HOUSE
            Matrix<T3> m(M+N,K);
            m.rowRange(0,N) = ma;
            m.rowRange(N,M+N) = ma;
            if (!(Norm(Hm-m) <= 0.001*Norm(m0)*Norm(H))) {
                cerr<<"Block2HouseholderLMult\n";
                cerr<<"Input: Y = "<<Y0<<endl;
                cerr<<"Z = "<<Z0<<endl;
                cerr<<"m = "<<m0<<endl;
                cerr<<"H = "<<H<<endl;
                cerr<<"Y = "<<Y<<endl;
                cerr<<"temp = "<<temp<<endl;
                cerr<<"Output: ma = "<<ma<<endl;
                cerr<<"mb = "<<mb<<endl;
                cerr<<"Hm = "<<Hm<<endl;
                abort();
            }
#endif
        }
    };

    // algo 90: Call inst
    template <class M1, class M2, class M3, class M4, class M5>
    struct Block2HouseholderLMult_Helper<90,M1,M2,M3,M4,M5>
    {
        static void call(const M1& Y, const M2& Z, M3& ma, M4& mb, M5& temp)
        {
            InstBlock2HouseholderLMult(
                Y.xView(),Z.xView(),ma.xView(),mb.xView(),temp.xView()); 
        }
    };

    // algo 95: Conjugate Z to match Y
    template <class M1, class M2, class M3, class M4, class M5>
    struct Block2HouseholderLMult_Helper<95,M1,M2,M3,M4,M5>
    {
        static void call(const M1& Y, const M2& Z, M3& ma, M4& mb, M5& temp)
        { 
            typedef typename M2::const_conjugate_type M2c;
            M2c Zc = Z.conjugate();
            Block2HouseholderLMult_Helper<-2,M1,M2c,M3,M4,M5>::call(
                Y,Zc,ma,mb,temp);
        }
    };

    // algo 97: Conjugate m
    template <class M1, class M2, class M3, class M4, class M5>
    struct Block2HouseholderLMult_Helper<97,M1,M2,M3,M4,M5>
    {
        static void call(const M1& Y, const M2& Z, M3& ma, M4& mb, M5& temp)
        { 
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            typedef typename M3::conjugate_type M3c;
            typedef typename M4::conjugate_type M4c;
            M1c Yc = Y.conjugate();
            M2c Zc = Z.conjugate();
            M3c mac = ma.conjugate();
            M4c mbc = mb.conjugate();
            Block2HouseholderLMult_Helper<-2,M1c,M2c,M3c,M4c,M5>::call(
                Yc,Zc,mac,mbc,temp);
        }
    };

    // algo -3: Select algorithm
    template <class M1, class M2, class M3, class M4, class M5>
    struct Block2HouseholderLMult_Helper<-3,M1,M2,M3,M4,M5>
    {
        static void call(const M1& Y, const M2& Z, M3& ma, M4& mb, M5& temp)
        { 
            const int algo = 11;
#ifdef PRINTALGO_HOUSE
            std::cout<<"Inline Block2HouseholderLMult\n";
            std::cout<<"Y = "<<Y<<std::endl;
            std::cout<<"Z = "<<Z<<std::endl;
            std::cout<<"ma = "<<ma<<std::endl;
            std::cout<<"mb = "<<mb<<std::endl;
            std::cout<<"temp = "<<temp<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            Block2HouseholderLMult_Helper<algo,M1,M2,M3,M4,M5>::call(
                Y,Z,ma,mb,temp);
#ifdef PRINTALGO_HOUSE
            std::cout<<"ma => "<<ma<<std::endl;
            std::cout<<"mb => "<<mb<<std::endl;
#endif
        }
    };

    // algo -2: Check for inst
    template <class M1, class M2, class M3, class M4, class M5>
    struct Block2HouseholderLMult_Helper<-2,M1,M2,M3,M4,M5>
    {
        static void call(const M1& Y, const M2& Z, M3& ma, M4& mb, M5& temp)
        { 
            typedef typename M3::real_type RT;
            const bool inst = 
                (M1::_colsize == TMV_UNKNOWN || M1::_colsize > 16) &&
                (M1::_rowsize == TMV_UNKNOWN || M1::_rowsize > 16) &&
                (M2::_size == TMV_UNKNOWN || M2::_size > 16) &&
                (M3::_colsize == TMV_UNKNOWN || M3::_colsize > 16) &&
                (M3::_rowsize == TMV_UNKNOWN || M3::_rowsize > 16) &&
                Traits<RT>::isinst;
            const int algo = 
                M4::_conj ? 97 : 
                M2::_conj != int(M1::_conj) ? 95 :
                inst ? 90 :
                -3;
            Block2HouseholderLMult_Helper<algo,M1,M2,M3,M4,M5>::call(
                Y,Z,ma,mb,temp);
        }
    };

    template <class M1, class M2, class M3, class M4, class M5>
    struct Block2HouseholderLMult_Helper<-1,M1,M2,M3,M4,M5>
    {
        static void call(const M1& Y, const M2& Z, M3& ma, M4& mb, M5& temp)
        {
            Block2HouseholderLMult_Helper<-2,M1,M2,M3,M4,M5>::call(
                Y,Z,ma,mb,temp); 
        }
    };

    template <class M1, class M2, class M3, class M4, class M5>
    static inline void InlineBlock2HouseholderLMult(
        const BaseMatrix_Rec<M1>& Y, const BaseMatrix_Tri<M2>& Z,
        BaseMatrix_Rec_Mutable<M3>& ma, BaseMatrix_Rec_Mutable<M4>& mb,
        BaseMatrix_Rec_Mutable<M5>& temp)
    {
        TMVStaticAssert((Sizes<M2::_size,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M3::_colsize,M2::_size>::same));
        TMVStaticAssert((Sizes<M4::_colsize,M1::_colsize>::same));
        TMVStaticAssert((Sizes<M3::_rowsize,M4::_rowsize>::same));
        TMVStaticAssert((Sizes<M5::_colsize,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M5::_rowsize,M3::_rowsize>::same));
        TMVAssert(Z.size() == Y.rowsize());
        TMVAssert(Y.rowsize() > 0);
        TMVAssert(Y.colsize() > 0);
        TMVAssert(ma.colsize() == Z.size());
        TMVAssert(mb.colsize() == Y.colsize());
        TMVAssert(ma.rowsize() == mb.rowsize());
        TMVAssert(temp.colsize() == Y.rowsize());
        TMVAssert(temp.rowsize() == ma.rowsize());
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        typedef typename M4::cview_type M4v;
        typedef typename M5::cview_type M5v;
        TMV_MAYBE_CREF(M1,M1v) Yv = Y.cView();
        TMV_MAYBE_CREF(M2,M2v) Zv = Z.cView();
        TMV_MAYBE_REF(M3,M3v) mav = ma.cView();
        TMV_MAYBE_REF(M4,M4v) mbv = mb.cView();
        TMV_MAYBE_REF(M5,M5v) tv = temp.cView();
        Block2HouseholderLMult_Helper<-3,M1v,M2v,M3v,M4v,M4v>::call(
            Yv,Zv,mav,mbv,tv);
    }

    template <class M1, class M2, class M3, class M4, class M5>
    static inline void Block2HouseholderLMult(
        const BaseMatrix_Rec<M1>& Y, const BaseMatrix_Tri<M2>& Z,
        BaseMatrix_Rec_Mutable<M3>& ma, BaseMatrix_Rec_Mutable<M4>& mb,
        BaseMatrix_Rec_Mutable<M5>& temp)
    {
        TMVStaticAssert((Sizes<M2::_size,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M3::_colsize,M2::_size>::same));
        TMVStaticAssert((Sizes<M4::_colsize,M1::_colsize>::same));
        TMVStaticAssert((Sizes<M3::_rowsize,M4::_rowsize>::same));
        TMVStaticAssert((Sizes<M5::_colsize,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M5::_rowsize,M3::_rowsize>::same));
        TMVAssert(Z.size() == Y.rowsize());
        TMVAssert(Y.rowsize() > 0);
        TMVAssert(Y.colsize() > 0);
        TMVAssert(ma.colsize() == Z.size());
        TMVAssert(mb.colsize() == Y.colsize());
        TMVAssert(ma.rowsize() == mb.rowsize());
        TMVAssert(temp.colsize() == Y.rowsize());
        TMVAssert(temp.rowsize() == ma.rowsize());
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        typedef typename M4::cview_type M4v;
        typedef typename M5::cview_type M5v;
        TMV_MAYBE_CREF(M1,M1v) Yv = Y.cView();
        TMV_MAYBE_CREF(M2,M2v) Zv = Z.cView();
        TMV_MAYBE_REF(M3,M3v) mav = ma.cView();
        TMV_MAYBE_REF(M4,M4v) mbv = mb.cView();
        TMV_MAYBE_REF(M5,M5v) tv = temp.cView();
        Block2HouseholderLMult_Helper<-2,M1v,M2v,M3v,M4v,M4v>::call(
            Yv,Zv,mav,mbv,tv);
    }


    //
    // Block2HouseholderLDiv
    //

    template <int algo, class M1, class M2, class M3, class M4, class M5>
    struct Block2HouseholderLDiv_Helper;

    // algo 11: Normal case
    template <class M1, class M2, class M3, class M4, class M5>
    struct Block2HouseholderLDiv_Helper<11,M1,M2,M3,M4,M5>
    {
        static void call(const M1& Y, const M2& Z, M3& ma, M4& mb, M5& temp)
        {
#ifdef XDEBUG_HOUSE
            int M = Y.colsize();
            int N = Y.rowsize();
            int K = ma.rowsize();

            typedef typename M1::value_type T;
            typedef typename M3::value_type T3;
            Matrix<T> Y0(M+N,N);
            Y0.rowRange(0,N).setToIdentity();
            Y0.rowRange(N,M+N) = Y;
            Matrix<T> Z0(Z);
            Matrix<T3> m0(M+N,K);
            m0.rowRange(0,N) = ma;
            m0.rowRange(N,N+M) = mb;
            Matrix<T> Hinv = T(1) - Y0*Z0.adjoint()*Y0.adjoint();
            Matrix<T3> Hm = Hinv*m0;
#endif

            // temp = -ZYtm
            // Remember that the full Y includes an identity matrix on top.
            temp = Y.adjoint() * mb;
            temp += ma;
            temp = -Z.adjoint() * temp;
            ma += temp;
            mb += Y * temp;

#ifdef XDEBUG_HOUSE
            Matrix<T3> m(M+N,K);
            m.rowRange(0,N) = ma;
            m.rowRange(N,M+N) = ma;
            if (!(Norm(Hm-m) <= 0.001*Norm(m0)*Norm(Hinv))) {
                cerr<<"Block2HouseholderLDiv\n";
                cerr<<"Input: Y = "<<Y0<<endl;
                cerr<<"Z = "<<Z0<<endl;
                cerr<<"m = "<<m0<<endl;
                cerr<<"Hinv = "<<Hinv<<endl;
                cerr<<"Y = "<<Y<<endl;
                cerr<<"temp = "<<temp<<endl;
                cerr<<"Output: ma = "<<ma<<endl;
                cerr<<"mb = "<<mb<<endl;
                cerr<<"Hm = "<<Hm<<endl;
                abort();
            }
#endif
        }
    };

    // algo 90: Call inst
    template <class M1, class M2, class M3, class M4, class M5>
    struct Block2HouseholderLDiv_Helper<90,M1,M2,M3,M4,M5>
    {
        static void call(const M1& Y, const M2& Z, M3& ma, M4& mb, M5& temp)
        {
            InstBlock2HouseholderLDiv(
                Y.xView(),Z.xView(),ma.xView(),mb.xView(),temp.xView()); 
        }
    };

    // algo 95: Conjugate Z to match Y
    template <class M1, class M2, class M3, class M4, class M5>
    struct Block2HouseholderLDiv_Helper<95,M1,M2,M3,M4,M5>
    {
        static void call(const M1& Y, const M2& Z, M3& ma, M4& mb, M5& temp)
        { 
            typedef typename M2::const_conjugate_type M2c;
            M2c Zc = Z.conjugate();
            Block2HouseholderLDiv_Helper<-2,M1,M2c,M3,M4,M5>::call(
                Y,Zc,ma,mb,temp);
        }
    };

    // algo 97: Conjugate m
    template <class M1, class M2, class M3, class M4, class M5>
    struct Block2HouseholderLDiv_Helper<97,M1,M2,M3,M4,M5>
    {
        static void call(const M1& Y, const M2& Z, M3& ma, M4& mb, M5& temp)
        { 
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            typedef typename M3::conjugate_type M3c;
            typedef typename M4::conjugate_type M4c;
            M1c Yc = Y.conjugate();
            M2c Zc = Z.conjugate();
            M3c mac = ma.conjugate();
            M4c mbc = mb.conjugate();
            Block2HouseholderLDiv_Helper<-2,M1c,M2c,M3c,M4c,M5>::call(
                Yc,Zc,mac,mbc,temp);
        }
    };

    // algo -3: Select algorithm
    template <class M1, class M2, class M3, class M4, class M5>
    struct Block2HouseholderLDiv_Helper<-3,M1,M2,M3,M4,M5>
    {
        static void call(const M1& Y, const M2& Z, M3& ma, M4& mb, M5& temp)
        { 
            const int algo = 11;
#ifdef PRINTALGO_HOUSE
            std::cout<<"Inline Block2HouseholderLDiv\n";
            std::cout<<"Y = "<<Y<<std::endl;
            std::cout<<"Z = "<<Z<<std::endl;
            std::cout<<"ma = "<<ma<<std::endl;
            std::cout<<"mb = "<<mb<<std::endl;
            std::cout<<"temp = "<<temp<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            Block2HouseholderLDiv_Helper<algo,M1,M2,M3,M4,M5>::call(
                Y,Z,ma,mb,temp);
#ifdef PRINTALGO_HOUSE
            std::cout<<"ma => "<<ma<<std::endl;
            std::cout<<"mb => "<<mb<<std::endl;
#endif
        }
    };

    // algo -2: Check for inst
    template <class M1, class M2, class M3, class M4, class M5>
    struct Block2HouseholderLDiv_Helper<-2,M1,M2,M3,M4,M5>
    {
        static void call(const M1& Y, const M2& Z, M3& ma, M4& mb, M5& temp)
        { 
            typedef typename M3::real_type RT;
            const bool inst = 
                (M1::_colsize == TMV_UNKNOWN || M1::_colsize > 16) &&
                (M1::_rowsize == TMV_UNKNOWN || M1::_rowsize > 16) &&
                (M2::_size == TMV_UNKNOWN || M2::_size > 16) &&
                (M3::_colsize == TMV_UNKNOWN || M3::_colsize > 16) &&
                (M3::_rowsize == TMV_UNKNOWN || M3::_rowsize > 16) &&
                Traits<RT>::isinst;
            const int algo = 
                M4::_conj ? 97 : 
                M2::_conj != int(M1::_conj) ? 95 :
                inst ? 90 :
                -3;
            Block2HouseholderLDiv_Helper<algo,M1,M2,M3,M4,M5>::call(
                Y,Z,ma,mb,temp);
        }
    };

    template <class M1, class M2, class M3, class M4, class M5>
    struct Block2HouseholderLDiv_Helper<-1,M1,M2,M3,M4,M5>
    {
        static void call(const M1& Y, const M2& Z, M3& ma, M4& mb, M5& temp)
        {
            Block2HouseholderLDiv_Helper<-2,M1,M2,M3,M4,M5>::call(
                Y,Z,ma,mb,temp); 
        }
    };

    template <class M1, class M2, class M3, class M4, class M5>
    static inline void InlineBlock2HouseholderLDiv(
        const BaseMatrix_Rec<M1>& Y, const BaseMatrix_Tri<M2>& Z,
        BaseMatrix_Rec_Mutable<M3>& ma, BaseMatrix_Rec_Mutable<M4>& mb,
        BaseMatrix_Rec_Mutable<M5>& temp)
    {
        TMVStaticAssert((Sizes<M2::_size,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M3::_colsize,M2::_size>::same));
        TMVStaticAssert((Sizes<M4::_colsize,M1::_colsize>::same));
        TMVStaticAssert((Sizes<M3::_rowsize,M4::_rowsize>::same));
        TMVStaticAssert((Sizes<M5::_colsize,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M5::_rowsize,M3::_rowsize>::same));
        TMVAssert(Z.size() == Y.rowsize());
        TMVAssert(Y.rowsize() > 0);
        TMVAssert(Y.colsize() > 0);
        TMVAssert(ma.colsize() == Z.size());
        TMVAssert(mb.colsize() == Y.colsize());
        TMVAssert(ma.rowsize() == mb.rowsize());
        TMVAssert(temp.colsize() == Y.rowsize());
        TMVAssert(temp.rowsize() == ma.rowsize());
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        typedef typename M4::cview_type M4v;
        typedef typename M5::cview_type M5v;
        TMV_MAYBE_CREF(M1,M1v) Yv = Y.cView();
        TMV_MAYBE_CREF(M2,M2v) Zv = Z.cView();
        TMV_MAYBE_REF(M3,M3v) mav = ma.cView();
        TMV_MAYBE_REF(M4,M4v) mbv = mb.cView();
        TMV_MAYBE_REF(M5,M5v) tv = temp.cView();
        Block2HouseholderLDiv_Helper<-3,M1v,M2v,M3v,M4v,M5v>::call(
            Yv,Zv,mav,mbv,tv);
    }

    template <class M1, class M2, class M3, class M4, class M5>
    static inline void Block2HouseholderLDiv(
        const BaseMatrix_Rec<M1>& Y, const BaseMatrix_Tri<M2>& Z,
        BaseMatrix_Rec_Mutable<M3>& ma, BaseMatrix_Rec_Mutable<M4>& mb,
        BaseMatrix_Rec_Mutable<M5>& temp)
    {
        TMVStaticAssert((Sizes<M2::_size,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M3::_colsize,M2::_size>::same));
        TMVStaticAssert((Sizes<M4::_colsize,M1::_colsize>::same));
        TMVStaticAssert((Sizes<M3::_rowsize,M4::_rowsize>::same));
        TMVStaticAssert((Sizes<M5::_colsize,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M5::_rowsize,M3::_rowsize>::same));
        TMVAssert(Z.size() == Y.rowsize());
        TMVAssert(Y.rowsize() > 0);
        TMVAssert(Y.colsize() > 0);
        TMVAssert(ma.colsize() == Z.size());
        TMVAssert(mb.colsize() == Y.colsize());
        TMVAssert(ma.rowsize() == mb.rowsize());
        TMVAssert(temp.colsize() == Y.rowsize());
        TMVAssert(temp.rowsize() == ma.rowsize());
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        typedef typename M4::cview_type M4v;
        typedef typename M5::cview_type M5v;
        TMV_MAYBE_CREF(M1,M1v) Yv = Y.cView();
        TMV_MAYBE_CREF(M2,M2v) Zv = Z.cView();
        TMV_MAYBE_REF(M3,M3v) mav = ma.cView();
        TMV_MAYBE_REF(M4,M4v) mbv = mb.cView();
        TMV_MAYBE_REF(M5,M5v) tv = temp.cView();
        Block2HouseholderLDiv_Helper<-2,M1v,M2v,M3v,M4v,M5v>::call(
            Yv,Zv,mav,mbv,tv);
    }

    // Finally a little helper function for determining the determinant
    // of the product of a bunch of Householder reflections.  When 
    // beta == 0, the determinant of H is 1, otherwise it is -1.
    // So from a list of beta's, we just *= -1 for each beta != 1.
    template <class V>
    int CalculateDetQ(const BaseVector<V>& beta)
    {
        typedef typename V::value_type RT;
        TMVStaticAssert(Traits<RT>::isreal);
        int detq = 1;
        const int n = beta.size();
        for(int i=0; i<n; ++i) if (beta[i] != RT(0)) detq = -detq;
        return detq;
    }


} // namespace tmv

#endif
