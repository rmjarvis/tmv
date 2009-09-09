//
// This file contains the code to implement Householder reflections.
//
// A Householder reflection is a unitary matrix H which transforms
// a given vector v into y e0 = (y,0,0,0...0) where |y| = |x|
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
// choosing a real y ( +-|x| ) or choosing y so that beta is real.
//
// The better choice is to make beta real, since then the matrix is
// Hermitian as well as unitary (ie. H = Ht, H*H = I), which also
// leads to the determinant being -1, which is actually required for
// it to be a "reflection" matrix.  It also has the advantage of (slightly)
// quicker multiplies, since the multiplication of H involves a multiplication
// by beta.  This step takes half the time if beta is real.
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
// However, I do want to take advantage of the LaPack optimized packages,
// so I bite the bullet and use their version of beta.  :P
//
// We still have two possible real values for y: y = +- |x|.
// We want to minimize the effects of rounding errors, which means maximize
// 1/beta = |x|^2 -+ |x| x0
// If real(x0) < 0, use y=-|x|, if real(x0) > 0, use y=|x|.
//
// The determinant of this Householder non-reflection is:
//
// det(H) = det( I - beta (x+|x|e0) (x+|x|e0)t )
//
// A simple, but not (to me) obvious identity is:
//
// det(I-beta v vt) = 1 - beta |v|^2
//
// Thus, det(H) = 1 - beta |x + |x|e0|^2 = -beta*/beta = -beta*^2/|beta|^2
//

#ifndef TMV_Householder_H
#define TMV_Householder_H

#include "TMV_Matrix.h"

namespace tmv {

  template <class T> T Householder_Reflect(const MatrixView<T>& M, T& det);
  // Multiplies M by a Householder matrix H which rotates
  // the first column into (y,0,0,0..0).
  // The vector v of the Householder matrix is stored
  // in the first column of M, except for the first element.
  // Tau is the return value. 
  // The rest of the matrix is multiplied by H.
  // (For the first column, this menas that y is the first element.)

  template <class T1, class T2> void Householder_Mult(
      const GenVector<T1>& v, T1 beta, const VectorView<T2>& z);
    // The input vector, v, is taken to be the vector for a  
    // Householder matrix, H.  This routine takes z <- H z
    // H z = z - beta v vt z

  template <class T1, class T2> void Householder_Mult(
      const GenVector<T1>& v, T1 beta, const MatrixView<T2>& M);
    // The input vector, v, is taken to be the vector for a  
    // Householder matrix, H.  This routine takes M <- H M
    // H M = M - beta v vt M

  template <class T1, class T2> void Householder_ConjMult(
      const GenVector<T1>& v, T1 beta, const VectorView<T2>& z)
  {
    // The input vector, v, is taken to be the vector for a  
    // Householder matrix, H.  This routine takes z <- H* z
    // H* z = z - beta* v* vT z
    TMVAssert(v.size() == z.size()-1);

    // Multiply z by H = I - beta* v* vT
    if (beta != T1(0)) Householder_Mult(v.Conjugate(),CONJ(beta),z);
  }

  template <class T1, class T2> void Householder_ConjMult(
      const GenVector<T1>& v, T1 beta, const MatrixView<T2>& M)
  {
    // The input vector, v, is taken to be the vector for a  
    // Householder matrix, H.  This routine takes M <- H* M
    // H* M = M - beta* v* vT M
    TMVAssert(v.size() == M.colsize()-1);

    // Multiply M by H = I - beta* v* vT
    if (beta != T1(0)) Householder_Mult(v.Conjugate(),CONJ(beta),M);
  }

  template <class T> void Householder_Unpack(const VectorView<T>& v, T beta);
    // The input vector, v, after the first element is taken to be the 
    // vector for a  Householder matrix, H. 
    //  v is set to Ht times (1,0,0...0)

  template <class T> void Householder_Unpack(const MatrixView<T>& m, T beta);
    // The input matrix is taken to have a Householder vector
    // stored in the first column.  This routine multiplies
    // the rest of the matrix by the Householder matrix Ht.
    // The First column is then set to Ht times (1,0,0...0)

} // namespace tmv

#endif
