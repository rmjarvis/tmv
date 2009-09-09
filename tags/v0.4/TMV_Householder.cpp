
#include "TMV.h"
#include "TMV_Householder.h"

namespace tmv {


  template <class T> T Householder_Reflect(
      const MatrixView<T>& m, T& det)
  {
    // Multiplies m by a Householder matrix H which rotates
    // the first column into (y,0,0,0..0).  
    // The vector v of the  Householder matrix is stored
    // in the first column of m, except for the first element.
    // Tau is the return value.
    // The rest of the matrix is multiplied by H.
    // (For the first column, this menas that y is the new first element.)
    TMVAssert(m.colsize() > 0);
    TMVAssert(m.rowsize() > 0);

    // Determine normx = |x|
    RealType(T) normsqx = NormSq(m.col(0,1,m.colsize()));
    const T x0 = m(0,0);

    // if all of x other than first element are 0, H is identity
    if (normsqx == RealType(T)(0) && IMAG(x0) == RealType(T)(0)) {
      return 0; 
      // Determinant in this case is 1 (H = I), so leave det alone.
    }

    // Finish calculating normx in usual case
    RealType(T) normsqx0 = NORM(x0);
    normsqx += normsqx0;
    RealType(T) normx = SQRT(normsqx);

    // beta = 1 / (|x|^2 + |x| x0)
    T beta = RealType(T)(1) / (normsqx + (REAL(x0) > 0 ? normx : -normx) * x0);

    // y = +- |x|
    T y =  REAL(x0) > 0 ? -normx : normx;

    // H = I - beta v vt
    // with v = x - y e0 in first column
    // Renormalize beta,v so that v(0) = 1
    T v0 = x0-y;
    RealType(T) normv0 = NORM(v0);
    beta *= normv0;
    T invv0 = RealType(T)(1)/v0;

    VectorView<T> v = m.col(0,1,m.colsize());
    v *= invv0;
    // Multiply the rest of m by H = I - beta v vt
    Householder_Mult(v,beta,m.Cols(1,m.rowsize()));

    m(0,0) = y;

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

  template <class T1, class T2> void Householder_Mult(
      const GenVector<T1>& v, T1 beta, const VectorView<T2>& z)
  {
    // The input vector, v, is taken to be the vector for a  
    // Householder matrix, H.  This routine takes z <- H z
    TMVAssert(v.size() == z.size()-1);

    // H z = z - beta v vt z
    if (beta != T1(0)) {

      T2 betavtz = z(0);
      VectorView<T2> zz = z.SubVector(1,z.size());
      betavtz += v.Conjugate()*zz;
      betavtz *= beta;

      z(0) -= betavtz;
      zz -= v*betavtz;
    }
  }

  template <class T1, class T2> void Householder_Mult(
      const GenVector<T1>& v, T1 beta, const MatrixView<T2>& m)
  {
    // The input vector, v, is taken to be the vector for a  
    // Householder matrix, H.  This routine takes m <- H m
    TMVAssert(m.colsize() > 0);
    TMVAssert(v.size() == m.colsize()-1);

    // H m = m - beta v vt m
    // We don't store the first element of v, which = 1, so:
    // H m = m - beta ( 1 ) ( 1 vt ) m
    //                ( v )
    // = m - beta ( 1 ) ( 1 vt ) ( m0 )
    //            ( v )          ( mx )
    // = m - beta ( 1 ) ( m0 + vt mx )
    //            ( v ) 
    if (m.rowsize() == 0) return;
    else if (m.rowsize() == 1) Householder_Mult(v,beta,m.col(0));
    else if (beta != T1(0)) {

      Vector<T2> vtm = m.row(0);
      MatrixView<T2> mx = m.Rows(1,m.colsize());
      vtm += v.Conjugate()*mx;

      m.row(0) -= beta * vtm;
      mx -= beta * v^vtm;
    }
  }
  
  template <class T> void Householder_Unpack(const VectorView<T>& v, T beta)
  {
    // The input vector, v, after the first element is taken to be the 
    // vector for a  Householder matrix, H. 
    //  v is set to Ht times (1,0,0...0)
    TMVAssert(v.size() > 0);

    if (beta == T(0)) {
      // then all but v(0) is already 0
      v(0) = T(1);
    } else {
      // v <- (I-beta* v vt) e0 = e0 - beta* v 
      v(0) = T(1)-CONJ(beta);
      v.SubVector(1,v.size()) *= -CONJ(beta);
    }
  }

  template <class T> void Householder_Unpack(const MatrixView<T>& m, T beta)
  {
    // The input matrix is taken to have a Householder vector
    // stored in the first column.  This routine multiplies
    // the rest of the matrix by the Householder matrix Ht.
    // The First column is then set to Ht times (1,0,0...0)
    TMVAssert(m.colsize() > 0);
    TMVAssert(m.rowsize() > 0);

    if (m.iscm()) {
      // Multiply the rest of m by Ht = I - beta* v vt
      Householder_Mult(m.col(0,1,m.colsize()),CONJ(beta),
	  m.SubMatrix(0,m.colsize(),1,m.rowsize()));
      // Finally, make the first column equal to H times e0
      Householder_Unpack(m.col(0),beta);
    } else {
      Householder_Mult(m.col(0,1,m.colsize()),CONJ(beta),
	  m.SubMatrix(0,m.colsize(),1,m.rowsize()));
      Householder_Unpack(m.col(0),beta);
    }
  }

#define InstFile "TMV_Householder.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


