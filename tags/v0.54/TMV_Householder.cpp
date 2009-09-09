
#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Householder.h"

namespace tmv {

  template <class T> T Householder_Reflect(T& x0, const VectorView<T>& x,
      T& det)
  {
    // Finds the Householder matrix H which rotates v into y e0.
    // The vector v of the Householder matrix is stored in v,
    // except for the first element.
    // Beta is the return value.

    // Determine normx = |x|
    RealType(T) normsqx = NormSq(x);

    // if all of x other than first element are 0, H is identity
    if (normsqx == RealType(T)(0) && IMAG(x0) == RealType(T)(0)) {
      return 0; 
      // Determinant in this case is 1 (H = I), so leave det alone.
    }

    // Finish calculating normx in usual case
    RealType(T) normsqx0 = NORM(x0);
    normsqx += normsqx0;
    RealType(T) normx = SQRT(normsqx);

    // y = +- |x|
    RealType(T) y =  REAL(x0) > 0 ? -normx : normx;

    // beta = 1 / (|x|^2 + |x| x0)
    T beta = RealType(T)(1) / (normsqx - y * x0);

    // H = I - beta v vt
    // with v = x - y e0 in first column
    // Renormalize beta,v so that v(0) = 1
    T v0 = x0-y;
    RealType(T) normv0 = NORM(v0);
    beta *= normv0;
    T invv0 = RealType(T)(1)/v0;

    x *= invv0;

    x0 = y;

    // Determinant of H = -beta^2/|beta|^2
    // But we are actually keeping track of the determinant of 
    // Q which is now multiplied by Ht.
    // The determinant of Ht = -conj(beta)^2/|beta|^2
    if (det != T(0)) {
      if (IMAG(beta) == RealType(T)(0))
	det = -det;
      else {
	det *= -CONJ(beta*beta)/NORM(beta);
      }
    }
    return beta;
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
    TMVAssert(IMAG(y) == RealType(T)(0));

    RealType(T) normsqx1 = NormSq(x);

    // if all of x other than first element are 0, H is identity
    if (normsqx1 == RealType(T)(0)) {
      beta = 0; 
      return true; 
      // Determinant in this case is 1 (H = I), so leave det alone.
    }

    RealType(T) normsqx = SQR(REAL(y));

    RealType(T) normsqx0 = normsqx - normsqx1;

    if (normsqx0 < RealType(T)(0)) {
      cout<<"y = "<<y<<" < Norm(x) = "<<SQRT(normsqx1)<<endl;
      cout<<"x = "<<x<<endl;
      return false;
    }
    // Same consideration on the +-: Maximize 1/beta.
    RealType(T) x0 =  REAL(y) > 0 ? -SQRT(normsqx0) : SQRT(normsqx0);

    // beta = 1 / (|x|^2 + |x| x0)
    beta = RealType(T)(1) / (normsqx - y * x0);

    // H = I - beta v vt
    // with v = x - y e0 in first column
    // Renormalize beta,v so that v(0) = 1
    T v0 = x0-y;
    RealType(T) normv0 = NORM(v0);
    beta *= normv0;
    T invv0 = RealType(T)(1)/v0;

    x *= invv0;

    y = x0;

    return true;
  }

  template <class T> T Householder_Reflect(const MatrixView<T>& m, T& det)
  {
    // Multiplies m by a Householder matrix H which rotates
    // the first column into y e0.
    // The vector v of the  Householder matrix is stored
    // in the first column of m, except for the first element.
    // Tau is the return value.
    // The rest of the matrix is multiplied by H.
    // (For the first column, this menas that y is the new first element.)
    TMVAssert(m.colsize() > 0);
    TMVAssert(m.rowsize() > 0);

    const VectorView<T> v = m.col(0,1,m.colsize());
    T beta;
    if (m.isconj()) {
      T m00 = m(0,0); 
      beta = Householder_Reflect(m00,v,det);
      m(0,0) = m00;
    } else {
      beta = Householder_Reflect(*m.ptr(),v,det);
    }
    if (beta != T(0)) Householder_LMult(v,beta,m.Cols(1,m.rowsize()));
    return beta;
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

    if (m0.size() == 0) return;
    else if (beta != T1(0)) {

      Vector<T2> temp = v.Conjugate() * mx;
      temp += m0;
      temp *= beta;

      m0 -= temp;
      mx -= v^temp;
    }
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
      T m00 = m(0,0);
      Householder_Unpack(m00,m.col(0,1,m.colsize()),beta);
      m(0,0) = m00;
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
    size_t M = Y.colsize();
    size_t N = Y.rowsize()-1; // # of cols already computed

    if (beta == T(0)) {
      Z.col(N,0,N+1).Zero();
    } else if (N == 0) {
      Z(0,0) = beta;
    } else {
      ConstVectorView<T> v = Y.col(N,N+1,M);
      VectorView<T> z = Z.col(N,0,N);
      z = Y.SubMatrix(N+1,M,0,N).Adjoint()*v;
      z += Y.row(N,0,N).Conjugate();
      z = -beta * Z.SubTriMatrix(0,N) * z;
      Z(N,N) = beta;
    }
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

    const size_t M = Y.colsize();
    const size_t N = Y.rowsize();

    if (N==1) {
      Z(0,0) = CONJ(beta(0));
    } else if (N==2) {
      Z(0,0) = CONJ(beta(0));
      Z(1,1) = CONJ(beta(1));
      T temp = Y.col(0,2,M).Conjugate()*Y.col(1,2,M);
      temp += CONJ(Y(1,0));
      Z(0,1) = -Z(0,0)*Z(1,1)*temp;
    } else {
      size_t j1 = (N+1)/2;
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

    size_t M = Y.colsize();
    size_t N = Y.rowsize();

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

    size_t M = Y.colsize();
    size_t N = Y.rowsize();

    if (m.iscm()) {
      Matrix<T2,ColMajor> ZtYtm = 
	LowerTriMatrixViewOf(Y.Rows(0,N),UnitDiag).Adjoint() * m.Rows(0,N);
      ZtYtm += Y.Rows(N,M).Adjoint() * m.Rows(N,M);
      ZtYtm = Z.Adjoint() * ZtYtm;
      m.Rows(0,N) -= LowerTriMatrixViewOf(Y.Rows(0,N),UnitDiag) * ZtYtm;
      m.Rows(N,M) -= Y.Rows(N,M) * ZtYtm;
    } else {
      Matrix<T2,RowMajor> ZtYtm = 
	LowerTriMatrixViewOf(Y.Rows(0,N),UnitDiag).Adjoint() * m.Rows(0,N);
      ZtYtm += Y.Rows(N,M).Adjoint() * m.Rows(N,M);
      ZtYtm = Z.Adjoint() * ZtYtm;
      m.Rows(0,N) -= LowerTriMatrixViewOf(Y.Rows(0,N),UnitDiag) * ZtYtm;
      m.Rows(N,M) -= Y.Rows(N,M) * ZtYtm;
    }
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

    size_t M = Y.colsize();
    size_t N = Y.rowsize();

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

#define InstFile "TMV_Householder.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


