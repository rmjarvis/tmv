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
// However, I do want the ability to call LaPack routines when available,
// so I bite the bullet and use their version of beta.  :P
//
// We still have two possible real values for y: y = +- |x|.
// We want to minimize the effects of rounding errors, which means maximize
// 1/beta = |x|^2 -+ |x| x0
// If real(x0) > 0, use y=-|x|, if real(x0) < 0, use y=|x|.
//
// The determinant of this Householder non-reflection is:
//
// det(H) = det( I - beta (x+|x|e0) (x+|x|e0)t )
//
// A simple, but not (to me) obvious identity is:
//
// det(I-alpha u ut) = 1 - alpha |u|^2
//
// Thus, det(H) = 1 - beta |x + |x|e0|^2 = -beta*/beta = -beta*^2/|beta|^2
//

#ifndef TMV_Householder_H
#define TMV_Householder_H

#include "tmv/TMV_Vector.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"

namespace tmv {

    template <typename T>
    T HouseholderReflect(T& x0, VectorView<T> x, T& det);
    // Finds the Householder matrix H which rotates (x0,x) into y e0.
    // The vector v of the Householder matrix is stored
    // in x, except for the first element, which is 1.
    // Beta is the return value.  y is returned as x0.

    template <typename T>
    inline T HouseholderReflect(ConjRef<T> x0, VectorView<T> x, T& det)
    {
        x0.getRef() = TMV_CONJ(x0.getRef());
        return HouseholderReflect(x0.getRef(),x,det);
        // x0r ends up real, so we don't need to care about the conjugation
        // of the output value.
    }

    template <typename T>
    inline T HouseholderReflect(
        VarConjRef<T> x0, VectorView<T> x, T& det)
    {
        if (x0.isconj()) x0.getRef() = TMV_CONJ(x0.getRef());
        return HouseholderReflect(x0.getRef(),x,det);
    }

    template <typename T>
    T HouseholderReflect(VectorView<T> x, T& det);

    template <typename T>
    T HouseholderReflect(MatrixView<T> m, T& det);
    // Multiplies m by a Householder matrix H which rotates
    // the first column into y e0.
    // The vector v of the  Householder matrix is stored
    // in the first column of m, except for the first element.
    // Beta is the return value.
    // The rest of the matrix is multiplied by H.
    // (For the first column, this menas that y is the new first element.)

    template <typename T>
    bool HouseholderUnReflect(T& y, VectorView<T> x, T& beta);
    // This is similar, except that the roles of y and x0 are swapped.
    // That is, the final rotated vector y e0 is presumed known, as is
    // the bulk of the input vector.
    // The unreflected value of the first element in the vector is returned
    // as y.
    // The rest of the vector is transformed into the Householder vector.
    // This is used for downdating QR decompositions.
    // The return value is true if successful, false if |y| < |x|.

    template <typename T>
    inline bool HouseholderUnReflect(
        ConjRef<T> x0, VectorView<T> x, T& beta)
    {
        TMVAssert(TMV_IMAG(x0.getRef()) == TMV_RealType(T)(0));
        return HouseholderUnReflect(x0.getRef(),x,beta);
    }

    template <typename T>
    inline bool HouseholderUnReflect(
        VarConjRef<T> x0, VectorView<T> x, T& beta)
    {
        TMVAssert(TMV_IMAG(x0.getRef()) == TMV_RealType(T)(0));
        return HouseholderUnReflect(x0.getRef(),x,beta);
    }

    template <typename T1, typename T2>
    void HouseholderLMult(
        const GenVector<T1>& v, T1 beta,
        VectorView<T2> m0, MatrixView<T2> mx);
    // The input vector, v, is taken to be the vector for a
    // Householder matrix, H.  This routine takes
    // ( m0 ) <- H ( m0 ) = [ ( m0 ) - beta ( 1 ) ( 1 vt ) ( m0 ) ]
    // ( mx )      ( mx )   [ ( mx )        ( v )          ( mx ) ]
    //
    // ( m0 ) -= beta (   m0 + vt mx   )
    // ( mx )         ( v (m0 + vt mx) )

    template <typename T1, typename T2>
    void HouseholderLMult(
        const GenVector<T1>& v, T1 beta, MatrixView<T2> m);

    template <typename T>
    void HouseholderUnpack(T& v0, VectorView<T> v, T beta);
    // The input vector, v, is taken to be the vector for a  Householder
    // matrix, H.
    // (v0,v) is set to Ht times e0.

    template <typename T>
    inline void HouseholderUnpack(
        ConjRef<T> v0, VectorView<T> v, T beta)
    {
        T vv = v0;
        HouseholderUnpack(vv,v,beta);
        v0 = vv;
    }

    template <typename T>
    inline void HouseholderUnpack(
        VarConjRef<T> v0, VectorView<T> v, T beta)
    {
        if (v0.isconj()) {
            T vv = v0;
            HouseholderUnpack(vv,v,beta);
            v0 = vv;
        } else {
            HouseholderUnpack(v0.getRef(),v,beta);
        }
    }

    template <typename T>
    void HouseholderUnpack(MatrixView<T> m, T beta);
    // The input matrix is taken to have a Householder vector
    // stored in the first column (not including the first element).
    // This routine multiplies the rest of the matrix by the Householder
    // matrix Ht.
    // The First column is then set to Ht times e0.

    template <typename T>
    void BlockHouseholderAugment(
        const GenMatrix<T>& Y, UpperTriMatrixView<T> Z, T beta);
    // All but the last columns of the input matrices, Y,Z are such that
    // I - Y'Z'Y't is a Block Householder matrix (the product of several
    // individual Householder matrices).
    // The last column of Y has the vector for the next Householder matrix.
    // Z is output such that:
    // I - YZYt = (I - Y'Z'Y't) (I - beta v vt)

    template <typename T>
    void BlockHouseholderMakeZ(
        const GenMatrix<T>& Y, UpperTriMatrixView<T> Z,
        const GenVector<T>& beta);
    // This routine calculates the Z component of the BlockHouseholder
    // formulation for Q.  Y contains the v's for the Householder matrices,
    // and beta contains the beta's.
    //
    // The output BlockHouseholder matrix I-YZYt is the product of
    // the Adjoints, H0t H1t ... HNt, since this is the product which
    // we actually use for most calculations.
    // If you want the product of the H's, input beta.Conjugate instead.
    //
    // Note that the Y matrix is really just the unit lower trapezoidal
    // component of the input Y.

    template <typename T, typename T2>
    void BlockHouseholderLMult(
        const GenMatrix<T>& Y, const GenUpperTriMatrix<T>& Z,
        MatrixView<T2> m);
    // The input Y,Z are such that (I - YZYt) is a block Householder matrix.
    // The upper square portion of Y is taken to be unit lower triangular.
    // ie. the diagonal and upper triangular portion are not referenced.
    // The routine then finds m <- (I - YZYt) m
    template <typename T, typename T2>
    void BlockHouseholderLDiv(
        const GenMatrix<T>& Y, const GenUpperTriMatrix<T>& Z,
        MatrixView<T2> m);
    // As above, but m <- (I - YZYt)^-1 m

    template <typename T>
    void BlockHouseholderUnpack(
        MatrixView<T> Y, const GenUpperTriMatrix<T>& Z,
        MatrixView<T> m);
    // This routine multiplies the rest of the matrix by the
    // BlockHouseholder matrix Ht, defined by Y,Z.
    // Then Y is unpacked in place.


} // namespace tmv

#endif
