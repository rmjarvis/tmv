///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2008                                                        //
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



#include "TMV_Householder.h"
#include "TMV_Vector.h"
#include "TMV_Matrix.h"
#include "TMV_TriMatrix.h"
#include "TMV_TriMatrixArith.h"
#include "TMV_MatrixArith.h"
#include "TMV_VectorArith.h"
#include "TMV_VIt.h"

//#define XDEBUG

#ifdef XDEBUG
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

#define RT RealType(T)

  template <class T> T Householder_Reflect(T& x0, const VectorView<T>& x,
      T& det)
  {
#ifdef XDEBUG
    //cout<<"Householder Reflect: x0 = "<<x0<<", x = "<<x<<endl;
    Vector<T> xx(x.size()+1);
    xx(0) = x0;
    xx.SubVector(1,xx.size()) = x;
#endif

    // Finds the Householder matrix H which rotates v into y e0.
    // The vector v of the Householder matrix is stored in v,
    // except for the first element.
    // Beta is the return value.

    // Determine normx = |x|
    RT normsqx = NormSq(x);
    //cout<<"normsqx = "<<normsqx<<endl;

    // if all of x other than first element are 0, H is identity
    if (normsqx == RT(0) && IMAG(x0) == RT(0)) {
      // Set them all to x explicitly in case underflow let to the 0.
      x.Zero();
      return 0; 
      // Determinant in this case is 1 (H = I), so leave det alone.
    }

    // Finish calculating normx in usual case
    RT normsqx0 = NORM(x0);
    normsqx += normsqx0;
    RT normx = SQRT(normsqx);
    //cout<<"normsqx0 = "<<normsqx0<<", normx = "<<normx<<endl;

    T beta;
    if (normsqx * Epsilon<T>() == RT(0)) {
      // Then we need to rescale, since underflow will cause rounding errors
      const RT eps = Epsilon<T>();
      // Epsilon is a pure power of 2, so no rounding errors from rescaling.
      x /= eps;
      x0 /= eps;
      beta = Householder_Reflect(x0,x,det);
      x0 *= eps;
    } else if (RT(1)/normsqx == RT(0)) {
      // Then we have overflow, so we need to rescale:
      const RT eps = Epsilon<T>();
      RT scale = eps; normx *= eps;
      while (normx > RT(1)) { scale *= eps; normx *= eps; }
      x *= scale;
      x0 *= scale;
      beta = Householder_Reflect(x0,x,det);
      x0 /= scale;
    } else {
      // The usual case with no rescaling required:

      // y = +- |x|
      RT y =  REAL(x0) > 0 ? -normx : normx;
      //cout<<"y = "<<y<<endl;

      // beta = 1 / (|x|^2 + |x| x0)
      // H = I - beta v vt
      // with v = x - y e0 in first column
      // Renormalize beta,v so that v(0) = 1
      T v0 = x0-y;
      RT normv0 = NORM(v0);
      beta = normv0 / (normsqx - y * x0);
      //cout<<"beta = "<<beta<<endl;
      T invv0 = RT(1)/v0;
      //cout<<"v0 = "<<v0<<", normv0 = "<<normv0<<", invv0 = "<<invv0<<endl;
      //cout<<"beta -> "<<beta<<endl;

      x *= invv0;
      //cout<<"x -> "<<x<<endl;

      x0 = y;
      //cout<<"x0 -> "<<x0<<endl;

      // Determinant of H = -beta^2/|beta|^2
      // But we are actually keeping track of the determinant of 
      // Q which is now multiplied by Ht.
      // The determinant of Ht = -conj(beta)^2/|beta|^2
      if (det != T(0)) {
	if (IMAG(beta) == RT(0))
	  det = -det;
	else 
	  det *= -CONJ(beta*beta)/NORM(beta);
      }

    }
#ifdef XDEBUG
    Vector<T> vv(xx.size());
    vv(0) = T(1);
    vv.SubVector(1,vv.size()) = x;
    Matrix<T> H = T(1)-beta*(vv^vv.Conjugate());
    // Check the following:
    // H * xx = (Norm(xx),0,0.0...)
    // x0 = Norm(xx)
    Vector<T> Hxx = H * xx;
    if (ABS(ABS(Hxx(0))-ABS(x0)) > 0.0001*ABS(x0) ||
	ABS(Norm(xx)-ABS(x0)) > 0.0001*ABS(x0) ||
	Norm(Hxx.SubVector(1,Hxx.size())) > 0.0001*ABS(x0)) {
      cerr<<"Householder Reflect:\n";
      cerr<<"Input: x = "<<xx<<endl;
      cerr<<"x0 = "<<x0<<endl;
      cerr<<"Norm(x) = "<<Norm(xx)<<endl;
      cerr<<"Output: v = "<<vv<<endl;
      cerr<<"beta = "<<beta<<endl;
      cerr<<"H = "<<H<<endl;
      cerr<<"xx = "<<xx<<endl;
      cerr<<"Hxx = "<<Hxx<<endl;
      cerr<<"abs(abs(hxx(0))-abs(x0)) = "<<ABS(ABS(Hxx(0))-ABS(x0))<<endl;
      cerr<<"abs(Norm(xx)-abs(x0)) = "<<ABS(Norm(xx))-ABS(x0)<<endl;
      cerr<<"Norm(Hxx(1,N)) = "<<Norm(Hxx.SubVector(1,Hxx.size()))<<endl;
      abort();
    }
#endif
    return beta;
  }

  template <class T> T Householder_Reflect(const MatrixView<T>& m, T& det)
  {
    // Multiplies m by a Householder matrix H which rotates
    // the first column into y e0.
    // The vector v of the  Householder matrix is stored
    // in the first column of m, except for the first element.
    // Beta is the return value.
    // The rest of the matrix is multiplied by H.
    // (For the first column, this means that y is the new first element.)
    TMVAssert(m.colsize() > 0);
    TMVAssert(m.rowsize() > 0);
#ifdef XDEBUG
    Matrix<T> m0(m);
#endif

    const VectorView<T> v = m.col(0,1,m.colsize());
    T beta;
    if (m.isconj()) {
      T m00 = CONJ(*m.cptr());
      beta = Householder_Reflect(m00,v,det);
#ifdef TMVFLDEBUG
      TMVAssert(m.ptr() >= m.first);
      TMVAssert(m.ptr() < m.last);
#endif
      *m.ptr() = CONJ(m00);
    } else {
      beta = Householder_Reflect(*m.ptr(),v,det);
    }
    if (beta != T(0)) Householder_LMult(v,beta,m.Cols(1,m.rowsize()));
#ifdef XDEBUG
    Vector<T> vv(m.colsize());
    vv(0) = T(1);
    vv.SubVector(1,vv.size()) = m.col(0,1,vv.size());
    Matrix<T> H = T(1)-beta*(vv^vv.Conjugate());
    //cout<<"H = "<<H<<endl;
    //cout<<"1-beta*vv^vvt = "<<(T(1)-beta*(vv^vv.Conjugate()))<<endl;
    Matrix<T> Hm = H * m0;
    Matrix<T> Hm2 = m;
    Hm2.col(0,1,vv.size()).Zero();
    if (Norm(Hm-Hm2) > 0.001*Norm(m0)) {
      cerr<<"Householder Reflect\n";
      cerr<<"Input: m0 = "<<m0<<endl;
      cerr<<"Output: vv = "<<vv<<endl;
      cerr<<"beta = "<<beta<<endl;
      cerr<<"m = "<<m<<endl;
      cerr<<"H = "<<H<<endl;
      cerr<<"vv^vvt = "<<(vv^vv.Conjugate())<<endl;
      cerr<<"beta*vv^vvt = "<<beta*(vv^vv.Conjugate())<<endl;
      cerr<<"1-beta*vv^vvt = "<<(T(1)-beta*(vv^vv.Conjugate()))<<endl;
      cerr<<"Hm = "<<Hm<<endl;
      cerr<<"Hm2 = "<<Hm2<<endl;
      abort();
    }
#endif
    return beta;
  }

  template <class T> T Householder_Reflect(ConjRef<T> x0,
      const VectorView<T>& x, T& det)
  {
    T& x0r = x0.GetRef();
    x0r = CONJ(x0r);
    return Householder_Reflect(x0r,x,det);
    // x0r ends up real, so we don't need to care about the conjugation
    // of the output value.
  }

  template <class T> T Householder_Reflect(VarConjRef<T> x0,
      const VectorView<T>& x, T& det)
  {
    T& x0r = x0.GetRef();
    if (x0.isconj()) x0r = CONJ(x0r);
    return Householder_Reflect(x0r,x,det);
  }

  template <class T> T Householder_Reflect(const VectorView<T>& x, T& det)
  {
    // Same as above, but takes (x0,x) to be contiguous
    return Householder_Reflect(x(0),x.SubVector(1,x.size()),det);
  }

  template <class T> bool Householder_UnReflect(T& y, const VectorView<T>& x,
      T& beta)
  {
    // This is similar, except that the roles of y and x0 are swapped.
    // That is, the final rotated vector y e0 is presumed known, as is 
    // the bulk of the input vector.  
    // The unreflected value of the first element in the vector is returned 
    // as y.
    // The rest of the vector is transformed into the Householder vector.
    // This is used for downdating QR decompositions.
    // The return value is true if successful, false if |y| < |x|.

    TMVAssert(x.size() > 0);
    TMVAssert(IMAG(y) == RT(0));

    RT normsqx1 = NormSq(x);

    // if all of x other than first element are 0, H is identity
    if (normsqx1 == RT(0)) {
      beta = 0; 
      return true; 
      // Determinant in this case is 1 (H = I), so leave det alone.
    }

    RT normsqx = SQR(REAL(y));
    RT normsqx0 = normsqx - normsqx1;
    if (normsqx0 < RT(0)) return false;

    // Same consideration on the +-: Maximize 1/beta.
    RT x0 =  REAL(y) > 0 ? -SQRT(normsqx0) : SQRT(normsqx0);

    // beta = 1 / (|x|^2 + |x| x0)
    // H = I - beta v vt
    // with v = x - y e0 in first column
    // Renormalize beta,v so that v(0) = 1
    T v0 = x0-y;
    RT normv0 = NORM(v0);
    beta = normv0 / (normsqx - y * x0);

    x /= v0;
    y = x0;

    return true;
  }

  template <class T> bool Householder_UnReflect(ConjRef<T> x0,
      const VectorView<T>& x, T& beta)
  {
    T& x0r = x0.GetRef();
    TMVAssert(IMAG(x0r) == RT(0));
    return Householder_UnReflect(x0r,x,beta);
  }

  template <class T> bool Householder_UnReflect(VarConjRef<T> x0,
      const VectorView<T>& x, T& beta)
  {
    TMVAssert(IMAG(x0.GetRef()) == RT(0));
    return Householder_UnReflect(x0.GetRef(),x,beta);
  }

  template <class T1, class T2> void Householder_LMult(
      const GenVector<T1>& v, T1 beta, const VectorView<T2>& m0,
      const MatrixView<T2>& mx)
  {
    // The input vector, v, is taken to be the vector for a  
    // Householder matrix, H.  This routine takes 
    // ( m0 ) <- H ( m0 ) = [ ( m0 ) - beta ( 1 ) ( 1 vt ) ( m0 ) ]
    // ( mx )      ( mx )   [ ( mx )        ( v )          ( mx ) ]
    // 
    // ( m0 ) -= beta (   m0 + vt mx   )
    // ( mx )         ( v (m0 + vt mx) )

    TMVAssert(v.size() == mx.colsize());
    TMVAssert(m0.size() == mx.rowsize());
#ifdef XDEBUG
    Matrix<T2> mm(mx.colsize()+1,m0.size());
    mm.row(0) = m0;
    mm.Rows(1,mm.colsize()) = mx;
#endif

    if (m0.size() == 0) return;
    else if (beta != T1(0)) {

      Vector<T2> temp = v.Conjugate() * mx;
      temp += m0;
      temp *= beta;

      m0 -= temp;
      mx -= v^temp;
    }
#ifdef XDEBUG
    Vector<T1> vv(v.size()+1);
    vv(0) = T1(1);
    vv.SubVector(1,vv.size()) = v;
    Matrix<T1> H = T1(1) - beta*(vv^vv.Conjugate());
    Matrix<T2> Hm = H * mm;
    Matrix<T2> Hm2(mx.colsize()+1,m0.size());
    Hm2.row(0) = m0;
    Hm2.Rows(1,Hm2.colsize()) = mx;
    if (Norm(Hm-Hm2) > 0.001*Norm(Hm)) {
      cerr<<"Householder_LMult\n";
      cerr<<"Input: m = "<<mm<<endl;
      cerr<<"vv = "<<vv<<endl;
      cerr<<"H = "<<H<<endl;
      cerr<<"Hm = "<<Hm<<endl;
      cerr<<"Output: m = "<<Hm2<<endl;
      abort();
    }
#endif
  }

  template <class T1, class T2> void Householder_LMult(
      const GenVector<T1>& v, T1 beta, const MatrixView<T2>& m)
  {
    // The same as above, except m0 and mx are a single contiguous
    // matrix, m.
    TMVAssert(m.colsize() > 0);
    Householder_LMult(v,beta,m.row(0),m.Rows(1,m.colsize()));
  }
  
  template <class T> void Householder_Unpack(T& v0,
      const VectorView<T>& v, T beta)
  {
    // The input matrix is taken to have a Householder vector
    // stored in the first column (not including the first element.   
    // This routine multiplies the rest of the matrix by the Householder   
    // matrix Ht.
    // The First column is then set to Ht times e0.

    if (beta == T(0)) {
      // then all but v(0) is already 0
      v0 = T(1);
    } else {
      // v <- (I-beta* (1 v)T (1 v*)) e0 = e0 - beta* v 
      v0 = T(1)-CONJ(beta);
      v *= -CONJ(beta);
    }
  }

  template <class T> static void Householder_Unpack(ConjRef<T> v0,
      const VectorView<T>& v, T beta)
  {
    T vv = v0;
    Householder_Unpack(vv,v,beta);
    v0 = vv;
  }

  template <class T> static void Householder_Unpack(VarConjRef<T> v0,
      const VectorView<T>& v, T beta)
  {
    if (v0.isconj()) {
      T vv = v0;
      Householder_Unpack(vv,v,beta);
      v0 = vv;
    } else {
      Householder_Unpack(v0.GetRef(),v,beta);
    }
  }


  template <class T> void Householder_Unpack(const MatrixView<T>& m, T beta)
  {
    // The input matrix is taken to have a Householder vector
    // stored in the first column (not including the first element.   
    // This routine multiplies the rest of the matrix by the Householder   
    // matrix Ht.
    // The First column is then set to Ht times e0.

    TMVAssert(m.colsize() > 0);
    TMVAssert(m.rowsize() > 0);

    // Multiply the rest of m by Ht = I - beta* v vt
    Householder_LMult(m.col(0,1,m.colsize()),CONJ(beta),
	m.SubMatrix(0,m.colsize(),1,m.rowsize()));
    // Finally, make the first column equal to Ht times e0
    if (m.isconj()) {
      T m00 = CONJ(*m.cptr());
      Householder_Unpack(m00,m.col(0,1,m.colsize()),beta);
#ifdef TMVFLDEBUG
      TMVAssert(m.ptr() >= m.first);
      TMVAssert(m.ptr() < m.last);
#endif
      *m.ptr() = CONJ(m00);
    } else {
      Householder_Unpack(*m.ptr(),m.col(0,1,m.colsize()),beta);
    }
  }

  template <class T> void BlockHouseholder_Augment(
      const GenMatrix<T>& Y, const UpperTriMatrixView<T>& Z, T beta)
  {
    // All but the last columns of the input matrices, Y,Z are such that
    // I - Y'Z'Y't is a Block Householder matrix (the product of several
    // individual Householder matrices).
    // The last column of Y has the vector for the next Householder matrix.
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

    TMVAssert(Z.size() == Y.rowsize());
    TMVAssert(Y.rowsize() > 0);
    TMVAssert(Y.colsize() > 0);
    TMVAssert(!Z.isconj());
    int M = Y.colsize();
    int N = Y.rowsize()-1; // # of cols already computed
#ifdef XDEBUG
    Matrix<T> Y0(Y);
    UpperTriMatrixViewOf(Y0).SetToIdentity();
    Matrix<T> Z0(Z);
    Matrix<T> H0 = T(1) - 
      Y0.Cols(0,N)*Z.SubTriMatrix(0,N)*Y0.Cols(0,N).Adjoint();
    Matrix<T> H1(Y.colsize(),Y.colsize());
    H1.SetToIdentity();
    H1.SubMatrix(N,M,N,M) -= beta * (Y0.col(N,N,M)^Y0.col(N,N,M).Conjugate());
    Matrix<T> H2 = H0*H1;
#endif

    if (beta == T(0)) {
      Z.col(N,0,N+1).Zero();
    } else if (N == 0) {
#ifdef TMVFLDEBUG
      TMVAssert(Z.ptr() >= Z.first);
      TMVAssert(Z.ptr() < Z.last);
#endif
      *Z.ptr() = beta;
    } else {
      ConstVectorView<T> v = Y.col(N,N+1,M);
      VectorView<T> z = Z.col(N,0,N);
      z = Y.SubMatrix(N+1,M,0,N).Adjoint()*v;
      z += Y.row(N,0,N).Conjugate();
      z = -beta * Z.SubTriMatrix(0,N) * z;
      // Z(N,N) = beta
#ifdef TMVFLDEBUG
      TMVAssert(Z.ptr()+N*(Z.stepi()+Z.stepj()) >= Z.first);
      TMVAssert(Z.ptr()+N*(Z.stepi()+Z.stepj()) < Z.last);
#endif
      *(Z.ptr() + N*(Z.stepi()+Z.stepj())) = beta;
    }
#ifdef XDEBUG
    Matrix<T> YY(Y);
    UpperTriMatrixViewOf(YY).SetToIdentity();
    Matrix<T> HH = T(1) - YY * Z * YY.Adjoint();
    if (Norm(HH-H2) > 0.001*Norm(H2)) {
      cerr<<"BlockHouseholder_Augment\n";
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

  template <class T> void BlockHouseholder_MakeZ(
      const GenMatrix<T>& Y, const UpperTriMatrixView<T>& Z, 
      const GenVector<T>& beta)
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
  {
    TMVAssert(Y.rowsize() == Z.rowsize());
    TMVAssert(Y.rowsize() == beta.size());
    TMVAssert(Y.colsize() >= Y.rowsize());
    TMVAssert(!Z.isconj());
    TMVAssert(!Y.isconj());
    TMVAssert(!beta.isconj());
#ifdef XDEBUG
    Matrix<T> Y0(Y);
    UpperTriMatrixViewOf(Y0).SetToIdentity();
    Matrix<T> Z0(Z);
    Vector<T> beta0 = beta;
    Matrix<T> Htot(Y.colsize(),Y.colsize());
    Htot.SetToIdentity();
    for(int i=0;i<int(Y.rowsize());i++) {
      Matrix<T> H(Y.colsize(),Y.colsize());
      H.SetToIdentity();
      H.SubMatrix(i,Y.colsize(),i,Y.colsize()) -= beta(i) * 
	(Y0.col(i,i,Y0.colsize()) ^ Y0.col(i,i,Y0.colsize()).Conjugate());
      Htot *= H.Adjoint();
    }
#endif

    const int M = Y.colsize();
    const int N = Y.rowsize();

    if (N==1) {
      // I - Y Z Yt = I - beta* v vt
      // Therefore:
      // Y.col(0) = v
      // Z(0,0) = beta*
#ifdef TMVFLDEBUG
      TMVAssert(Z.ptr() >= Z.first);
      TMVAssert(Z.ptr() < Z.last);
#endif
      *Z.ptr() = CONJ(*beta.cptr());
    } else if (N==2) {
      // Y = ( Y00  0  )   Z = ( Z00  Z01 )
      //     ( Yx0 Yx1 )       (  0   Z11 )
      //
      // I - Y Z Yt = I - ( Z00 Y00   Z01 Y00           ) ( Y00* Yx0t )
      //                  ( Z00 Yx0   Z01 Yx0 + Z11 Yx1 ) (  0   Yx1t )
      //   = ( 1 - Z00 Y00 Y00*   -Z00 Y00 Yx0t - Z01 Y00 Yx1t               )
      //     ( -Z00 Y00* Yx0      I - Z00 Yx0 Yx0t - Z01 Yx0 Yx1t 
      //                                 - Z11 Yx1 Yx1t                      )
      //
      //
      // (I - b0* v0 v0t)*( 1            0     )
      //                  ( 0   I - b1* v1 v1t )
      // [ Let v0 = ( v00 v0x ) ]
      //    = ( 1 - b0* v00 v00*   -b0* v00 v0xt    ) ( 1        0        )
      //      ( -b0* v00* v0x      I - b0* v0x v0xt ) ( 0  I - b1* v1 v1t )
      //    = ( 1 - b0* v00 v00*   -b0* v00 v0xt + b0* b1* v00 v0xt v1 v1t )
      //      ( -b0* v00* v0x      I - b0* v0x v0xt - b1* v1 v1t 
      //                               + b0* b1* v0x v0xt v1 v1t           )
      //
      // Matching the two results, we find:
      //
      // Y00 = v00
      // Yx0 = v0x
      // Yx1 = v1
      // Z00 = b0*
      // Z11 = b1*
      // Z01 = -b0* b1* v0xt v1
      T* Z00 = Z.ptr();
      T* Z01 = Z00 + Z.stepj();
      T* Z11 = Z01 + Z.stepi();
      const T cb0 = CONJ(*beta.cptr());
      const T cb1 = CONJ(*(beta.cptr() + beta.step()));
      const T cY10 = CONJ(*(Y.cptr() + Y.stepi()));
      // Z(0,0) = CONJ(beta(0))
#ifdef TMVFLDEBUG
      TMVAssert(Z00 >= Z.first);
      TMVAssert(Z00 < Z.last);
      TMVAssert(Z01 >= Z.first);
      TMVAssert(Z01 < Z.last);
      TMVAssert(Z11 >= Z.first);
      TMVAssert(Z11 < Z.last);
#endif
      *Z00 = cb0;
      // Z(1,1) = CONJ(beta(1))
      *Z11 = cb1;
      T temp = Y.col(0,2,M).Conjugate()*Y.col(1,2,M);
      // temp += CONJ(Y(1,0))
      temp += cY10;
      // Z(0,1) = -Z(0,0)*Z(1,1)*temp;
      *Z01 = -cb0*cb1*temp;
    } else {
      int j1 = (N+1)/2;
      ConstMatrixView<T> Y1 = Y.Cols(0,j1);
      UpperTriMatrixView<T> Z1 = Z.SubTriMatrix(0,j1);
      ConstVectorView<T> beta1 = beta.SubVector(0,j1);
      BlockHouseholder_MakeZ(Y1,Z1,beta1);

      ConstMatrixView<T> Y2 = Y.SubMatrix(j1,M,j1,N);
      UpperTriMatrixView<T> Z2 = Z.SubTriMatrix(j1,N);
      ConstVectorView<T> beta2 = beta.SubVector(j1,N);
      BlockHouseholder_MakeZ(Y2,Z2,beta2);

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
      MatrixView<T> Z3 = Z.SubMatrix(0,j1,j1,N);
      Z3 = Y1.Rows(j1,N).Adjoint() *
	LowerTriMatrixViewOf(Y.SubMatrix(j1,N,j1,N),UnitDiag);
      Z3 += Y1.Rows(N,M).Adjoint() * Y.SubMatrix(N,M,j1,N);
      Z3 = -Z1*Z3;
      Z3 *= Z2;
    }
#ifdef XDEBUG
    Matrix<T> Y2(Y);
    UpperTriMatrixViewOf(Y2).SetToIdentity();
    Matrix<T> Hnet = T(1) - Y2*Z*Y2.Adjoint();
    if (Norm(Htot-Hnet) > 0.001*Norm(Htot)) {
      cerr<<"BlockHouseholder_MakeZ\n";
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

  template <class T, class T2> void BlockHouseholder_LMult(
      const GenMatrix<T>& Y, const GenUpperTriMatrix<T>& Z,
      const MatrixView<T2>& m)
  {
    // The input Y,Z are such that (I - YZYt) is a block Householder matrix.
    // The upper square portion of Y is taken to be unit lower triangular.
    // ie. the diagonal and upper triangular portion are not referenced.
    // The routine then finds m <- (I - YZYt) m

    TMVAssert(Z.size() == Y.rowsize());
    TMVAssert(Y.rowsize() > 0);
    TMVAssert(Y.colsize() > 0);
    TMVAssert(m.colsize() == Y.colsize());
#ifdef XDEBUG
    Matrix<T> Y0(Y);
    UpperTriMatrixViewOf(Y0).SetToIdentity();
    Matrix<T> Z0(Z);
    Matrix<T2> m0(m);
    Matrix<T> H = T(1) - Y0*Z0*Y0.Adjoint();
    Matrix<T2> Hm = H*m0;
#endif

    int M = Y.colsize();
    int N = Y.rowsize();

    if (m.iscm()) {
      Matrix<T2,ColMajor> ZYtm = 
	LowerTriMatrixViewOf(Y.Rows(0,N),UnitDiag).Adjoint() * m.Rows(0,N);
      ZYtm += Y.Rows(N,M).Adjoint() * m.Rows(N,M);
      ZYtm = Z * ZYtm;
      m.Rows(0,N) -= LowerTriMatrixViewOf(Y.Rows(0,N),UnitDiag) * ZYtm;
      m.Rows(N,M) -= Y.Rows(N,M) * ZYtm;
    } else {
      Matrix<T2,RowMajor> ZYtm = 
	LowerTriMatrixViewOf(Y.Rows(0,N),UnitDiag).Adjoint() * m.Rows(0,N);
      ZYtm += Y.Rows(N,M).Adjoint() * m.Rows(N,M);
      ZYtm = Z * ZYtm;
      m.Rows(0,N) -= LowerTriMatrixViewOf(Y.Rows(0,N),UnitDiag) * ZYtm;
      m.Rows(N,M) -= Y.Rows(N,M) * ZYtm;
    }
#ifdef XDEBUG
    if (Norm(Hm-m) > 0.001*Norm(m0)*Norm(H)) {
      cerr<<"BlockHouseholder_LMult\n";
      cerr<<"Input: Y = "<<Y0<<endl;
      cerr<<"Z = "<<Z0<<endl;
      cerr<<"m = "<<m0<<endl;
      cerr<<"H = "<<H<<endl;
      cerr<<"Output: m = "<<m<<endl;
      abort();
    }
#endif
  }

  template <class T, class T2> void BlockHouseholder_LDiv(
      const GenMatrix<T>& Y, const GenUpperTriMatrix<T>& Z,
      const MatrixView<T2>& m)
  {
    // The routine finds m <- (I - YZYt)^-1 m
    // = (I - YZtYt) m

    TMVAssert(Z.size() == Y.rowsize());
    TMVAssert(Y.rowsize() > 0);
    TMVAssert(Y.colsize() > 0);
    TMVAssert(m.colsize() == Y.colsize());
#ifdef XDEBUG
    Matrix<T> Y0(Y);
    UpperTriMatrixViewOf(Y0).SetToIdentity();
    Matrix<T> Z0(Z.Adjoint());
    Matrix<T2> m0(m);
    Matrix<T> Hinv = T(1) - Y0*Z0*Y0.Adjoint();
    Matrix<T2> Hm = Hinv*m0;
#endif

    int M = Y.colsize();
    int N = Y.rowsize();

    if (m.isrm()) {
      Matrix<T2,RowMajor> ZtYtm = 
	LowerTriMatrixViewOf(Y.Rows(0,N),UnitDiag).Adjoint() * m.Rows(0,N);
      ZtYtm += Y.Rows(N,M).Adjoint() * m.Rows(N,M);
      ZtYtm = Z.Adjoint() * ZtYtm;
      m.Rows(0,N) -= LowerTriMatrixViewOf(Y.Rows(0,N),UnitDiag) * ZtYtm;
      m.Rows(N,M) -= Y.Rows(N,M) * ZtYtm;
    } else {
      Matrix<T2,ColMajor> ZtYtm = 
	LowerTriMatrixViewOf(Y.Rows(0,N),UnitDiag).Adjoint() * m.Rows(0,N);
      ZtYtm += Y.Rows(N,M).Adjoint() * m.Rows(N,M);
      ZtYtm = Z.Adjoint() * ZtYtm;
      m.Rows(0,N) -= LowerTriMatrixViewOf(Y.Rows(0,N),UnitDiag) * ZtYtm;
      m.Rows(N,M) -= Y.Rows(N,M) * ZtYtm;
    }
#ifdef XDEBUG
    if (Norm(Hm-m) > 0.001*Norm(m0)*Norm(Hinv)) {
      cerr<<"BlockHouseholder_LDiv\n";
      cerr<<"Input: Y = "<<Y0<<endl;
      cerr<<"Z = "<<Z0<<endl;
      cerr<<"m = "<<m0<<endl;
      cerr<<"Hinv = "<<Hinv<<endl;
      cerr<<"Output: m = "<<m<<endl;
      abort();
    }
#endif
  }

  template <class T> void BlockHouseholder_Unpack(
      const MatrixView<T>& Y, const GenUpperTriMatrix<T>& Z,
      const MatrixView<T>& m)
  {
    // This routine multiplies the rest of the matrix by the 
    // BlockHouseholder matrix Ht, defined by Y,Z.
    // Then Y is unpacked in place.
    TMVAssert(Y.colsize() > 0);
    TMVAssert(Y.rowsize() > 0);
    TMVAssert(Y.colsize() >= Y.rowsize());
    TMVAssert(Y.colsize() == m.colsize());
    TMVAssert(Y.rowsize() == Z.size());

    int M = Y.colsize();
    int N = Y.rowsize();

    // Multiply the rest of m by Ht
    BlockHouseholder_LMult(Y,Z,m);
    // Make the first N columns equal to 
    // Ht [ I ] = (I - YZYt) [ I ]
    //    [ 0 ]              [ 0 ]
    UpperTriMatrix<T,NonUnitDiag,RowMajor> temp = 
      -Z * LowerTriMatrixViewOf(Y.Rows(0,N),UnitDiag).Adjoint();
    Y.Rows(N,M) *= temp;
    Y.Rows(0,N) = LowerTriMatrixViewOf(Y.Rows(0,N),UnitDiag) * temp;
    Y.Rows(0,N).diag().AddToAll(T(1));
  }

#undef RT

#define InstFile "TMV_Householder.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


