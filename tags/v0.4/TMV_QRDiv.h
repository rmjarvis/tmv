//---------------------------------------------------------------------------
//
// This file contains the code for doing division using 
// QR Decomposition.
//
// The basic idea of an QR decomposition is that any 
// matrix A can be decomposed into a unitary matrix
// times an upper triangle matrix.
//
// A = Q R
//
// We do this using Householder transformations, which can
// be stored in a lower triangle matrix, thus other than the 
// diagonal, which they both need, we can store them in the 
// place of the original matrix. It is more convenient to 
// keep the diagonal of Q in place and take the diagonal of R
// separately.
//
// If R is not singular, the solution to A x = b is found from
// Q R x = b
// R x = Qt b
// which can be solved by back-substitutaion.
//
// If m > n, this does not actually give a solution to A x = b,
// since Q Qt != I  (Q is only column orthogonal if m > n.)
// But it does give the value of x which minimizes the 2-norm
// of (A x - b)
//
// The 2-norm of v is the square root of vt v, so
// |A x - b|^2 =
// (A x - b)t (A x - b) = 
// (xt At - bt) (A x - b) =
// xt At A x - bt A x - xt At b + bt b =
// xt Rt Qt Q R x - bt Q R x - xt Rt Qt b + bt b =
// |R x - Qt b|^2 + |b|^2 - |Qt b|^2
// Clearly the x which minimizes this is the solution of R x = Qt b.
//
// If R is singular, then back-substitution will fail in solving
// R x = Qt b.  However, there is a trick that still gives us a valid 
// least squares solution to A x = b.  This is to use column
// pivoting to obtain a decomposition:
// A = Q [ R11 R12 ] P
//       [  0   0  ]
// where R11 is upper triangular. 
//
// With this decomposition, 
// Q R P x = b
// R P x = Qt b
// Let z = P x (ie. we will solve for z first, then x = z / P = P * z)
// and c = Qt b
// Then R z = c
// Since R is singular, there is no exact solution, but for the least-squares
// problem, we just want to minimize |R z - c|.
// If z = [ z1 ] and c = [ c1 ] (with the same dimensional split as R11 and R12)
//        [ z2 ]         [ c2 ]
// then R z - c = [ R11 z1 + R12 z2 - c1 ]
//                [        -c2           ]
// The minimum norm will occur for any z2 if z1 = R11^-1 (c1 - R12 z2)
// So a solution can be found with z2 = 0.
// This solution then has minimum |A x - b| (ie. the least squares
// solution).  
//
// For dividing from the other side with m > n, we can always find a 
// possible solution (assuming R is non-singular):
//
// xt Q R P = bt 
// xt Q R = bt Pt
// xt Q = bt Pt R^-1  (This is done by front-substitution)
// xt = bt Pt R^-1 Qt
//
// This solution does satisfy the equation xt A = bt, but
// it will not be the solution with minimum |x|.  
// Use SVD for the min |x| solution.
//
// If R is singular, we of course have a problem with the front-substitution
// step.  This time, the QRP decomposition does not work quite as well.
// Take the front-substitution step (the Q and P steps don't affect 
// this minimization), and write R as above:
//
// zt R = ct
// [ z1t z2t ] [ R11  R12 ] = [ c1t c2t ]
//             [  0    0  ]
// [ (z1t R11) (z1t R12) ] = [ c1t c2t]
// |zt R - ct|^2 = |z1t R11 - c1t|^2 + |z1t R12 - c2t|^2
// We can set z2t = 0, since it is arbitrary.
// But it is not so clear what the correct z1t is in this case.
// We take z1t = c1t R11^-1, as an approximate solution, but point out 
// that this is not really correct. You should use SVD instead for 
// the correct least squares solution in this case.
//
// We provide two Divider classes here, QRDiv and QRPDiv.
// QRDiv does not include all of the column permutations, so it is
// a big faster, but it will fail for singular matrices.
//


#ifndef TMV_QRDiv_H
#define TMV_QRDiv_H

#include "TMV_Vector.h"
#include "TMV_Matrix.h"
#include "TMV_Permutation.h"
#include "TMV_PermutationArith.h"
#include "TMV_Divider.h"

namespace tmv {

  // These are defined in TMV_QRDiv.cpp
  
  template <class T> void QR_Decompose(
      const MatrixView<T>& QRx, const VectorView<T>& Qbeta, T& det);
  template <class T> void QR_Decompose(
      const MatrixView<T>& QRx, const VectorView<T>& Qbeta)
  { T d=0; QR_Decompose(QRx,Qbeta,d); }
  // Decompose A (input as QRx) into Q R.
  // On output, Q is stored in the lower triangle part of QRx as
  // Householder vectors, and the Qbeta vector.
  // R is the upper triangle part of QRx
 
  template <class T> void QRP_Decompose(
      const MatrixView<T>& QRx,
      const VectorView<T>& Qbeta, Permutation& P, T& det);
  template <class T> void QRP_Decompose(
      const MatrixView<T>& QRx,
      const VectorView<T>& Qbeta, Permutation& P)
  { T d=0; QRP_Decompose(QRx,Qbeta,P,d); }
  // Decompose A (input as QRx) into Q R P.
  // On output, Q is stored in the lower triangle part of QRx as
  // Householder vectors.  The upper triangle part and Qbeta hold R.
 
  template <class T> void QR_Decompose(
      const MatrixView<T>& Q, const MatrixView<T>& R, T& det);
  template <class T> void QR_Decompose(
      const MatrixView<T>& Q, const MatrixView<T>& R)
  { T d=0; QR_Decompose(Q,R,d); }
  // Decompose A (input as Q) into Q R.
 
  template <class T> void QRP_Decompose(
      const MatrixView<T>& Q, const MatrixView<T>& R, 
      Permutation& P, T& det);
  template <class T> void QRP_Decompose(
      const MatrixView<T>& Q, 
      const MatrixView<T>& R, Permutation& P)
  { T d=0; QRP_Decompose(Q,R,P,d); }
  // Decompose A (input as Q) into Q R P.
 
  template <class T> void GetQFromQR(const MatrixView<T>& Q,
      const GenVector<T>& Qbeta);
  template <class T1, class T2> void Q_LDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& Qbeta,
      const MatrixView<T2>& m);
  template <class T1, class T2> void Q_RDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& Qbeta,
      const MatrixView<T2>& m);

  template <class T> bool QR_DownDate(
      const MatrixView<T>& R, const GenVector<T>& z);
  // Given A = QR, where one of the rows of A is z,
  // find the QR decomposition for A' which is A with the row z removed.
  // This function only updates R, ignoring the Q matrix.

  template <class T1, class T2, class T3> void QR_LDiv(
      const GenMatrix<T1>& QR, const GenVector<T1>& Qbeta,
      const GenMatrix<T2>& m, const MatrixView<T3>& x, size_t N1=0);
  template <class T1, class T2> void QR_LDivEq(
      const GenMatrix<T1>& QR, const GenVector<T1>& Qbeta, 
      const MatrixView<T2>& m, size_t N1=0);
  template <class T1, class T2, class T3> void QR_RDiv(
	const GenMatrix<T1>& QR, const GenVector<T1>& Qbeta, 
	const GenMatrix<T2>& m, const MatrixView<T3>& x, size_t N1=0);
  template <class T1, class T2> void QR_RDivEq(
      const GenMatrix<T1>& QR, const GenVector<T1>& Qbeta, 
      const MatrixView<T2>& m, size_t N1=0);

#define CT complex<T>
  template <class T> void QR_LDivEq(
      const GenMatrix<CT>& QR, const GenVector<CT>& Qbeta, 
      const MatrixView<T>& m, size_t N1=0)
  { TMVAssert(false); }
  template <class T> void QR_RDivEq(
      const GenMatrix<CT>& QR, const GenVector<CT>& Qbeta, 
      const MatrixView<T>& m, size_t N1=0)
  { TMVAssert(false); }
#undef CT

  template <class T> class QRDiv : 
    virtual public BaseQRDiv<T> 
  {

    public :

      //
      // Constructors
      //

#define AA (istrans?A.QuickTranspose():A.QuickView())
      QRDiv(const GenMatrix<T>& A) : 
	istrans(A.isrm()),
	isshort(A.colsize()<A.rowsize()),
	QRx(AA), Qbeta(min(A.rowsize(),A.colsize())), 
	det(T(1)), donedet(false), N1(Qbeta.size())
      { 
	if (isshort!=istrans)
	  QR_Decompose(QRx.QuickTranspose(),Qbeta.View(),det); 
	else
	  QR_Decompose(QRx.QuickView(),Qbeta.View(),det); 
	while(N1>0 && QRx.diag()(N1-1)==T(0)) --N1;
      }
#undef AA

      ~QRDiv() {}

      //
      // Div, DivEq
      //

      template <class T1> void DoLDivEq(const MatrixView<T1>& m) const
      {
	TMVAssert(QRx.IsSquare());
	TMVAssert(m.colsize() == QRx.colsize());
	TMVAssert(!isshort);
	if (istrans)
	  QR_LDivEq(QRx.QuickTranspose(),Qbeta,m,N1);
	else
	  QR_LDivEq(QRx,Qbeta,m,N1);
      } 

      template <class T1> void DoRDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(QRx.IsSquare());
	TMVAssert(m.rowsize() == QRx.colsize());
	TMVAssert(!isshort);
	if (istrans)
	  QR_RDivEq(QRx.QuickTranspose(),Qbeta,m,N1);
	else
	  QR_RDivEq(QRx,Qbeta,m,N1);
      } 

      template <class T1, class T2> void DoLDiv(const GenMatrix<T1>& m, 
	  const MatrixView<T2>& x) const 
      {
	TMVAssert(m.rowsize() == x.rowsize());
	TMVAssert(m.colsize() == 
	    (isshort==istrans) ? QRx.colsize() : QRx.rowsize());
	TMVAssert(x.colsize() == 
	    (isshort==istrans) ? QRx.rowsize() : QRx.colsize());
	if (isshort) 
	  if (istrans)
	    QR_RDiv(QRx,Qbeta,m.QuickTranspose(),x.QuickTranspose(),N1);
	  else
	    QR_RDiv(QRx.QuickTranspose(),Qbeta,m.QuickTranspose(),
		x.QuickTranspose(),N1);
	else 
	  if (istrans)
	    QR_LDiv(QRx.QuickTranspose(),Qbeta,m,x,N1);
	  else
	    QR_LDiv(QRx,Qbeta,m,x,N1);
      } 

      template <class T1, class T2> void DoRDiv(const GenMatrix<T1>& m, 
	  const MatrixView<T2>& x) const 
      {
	TMVAssert(m.colsize() == x.colsize());
	TMVAssert(m.colsize() == 
	    (isshort==istrans) ? QRx.rowsize() : QRx.colsize());
	TMVAssert(x.colsize() == 
	    (isshort==istrans) ? QRx.colsize() : QRx.rowsize());
	if (isshort)
	  if (istrans) 
	    QR_LDiv(QRx,Qbeta,m.QuickTranspose(),x.QuickTranspose(),N1);
	  else
	    QR_LDiv(QRx.QuickTranspose(),Qbeta,m.QuickTranspose(),
		x.QuickTranspose(),N1);
	else 
	  if (istrans) 
	    QR_RDiv(QRx.QuickTranspose(),Qbeta,m,x,N1);
	  else
	    QR_RDiv(QRx,Qbeta,m,x,N1);
      } 

#include "TMV_AuxAllDiv.h"

      //
      // Determinant, Inverse
      //

      T Det() const;
      Matrix<T,ColMajor> Inverse() const;
      Matrix<T,ColMajor> DoInverseATA() const;
      inline bool Singular() const { return Det() == T(0); }
      inline Matrix<T,ColMajor> InverseATA() const
      {
#ifdef TMVDEBUG
	if (isshort) {
	  cout<<"Warning InverseATA called for short matrix in QRDiv\n";
	  cout<<"The result is really (AAt)^-1\n";
	}
#endif
	return DoInverseATA();
      }


      //
      // Access Decomposition
      //

      Matrix<T,ColMajor> QR_GetQ() const;
      Matrix<T,ColMajor> QR_GetR() const;
      inline Permutation QR_GetP() const { return Permutation(Qbeta.size()); }
      inline bool QR_IsTrans() const { return isshort; }
      const Matrix<T,ColMajor>& GetQR() const { return QRx; }
      const Vector<T>& GetQbeta() const { return Qbeta; }

      inline std::string Type() const
      { return std::string("QRDiv<") + tmv::Type(T()) + ">"; }

    protected :
      bool istrans;
      bool isshort;
      Matrix<T,ColMajor> QRx;
      Vector<T> Qbeta;
      mutable T det;
      mutable bool donedet;
      size_t N1;
  };

  template <class T> class QRPDiv : 
    virtual public BaseQRDiv<T> 
  {

    public :

      //
      // Constructors
      //

#define AA (istrans ? A.QuickTranspose() : A.QuickView())
      QRPDiv(const GenMatrix<T>& A) : 
	istrans(A.colsize() < A.rowsize()),
	QRx(AA), Qbeta(min(QRx.colsize(),QRx.rowsize())), 
	P(Qbeta.size()), det(T(1)), donedet(false), N1(Qbeta.size())
      { 
	QRP_Decompose(QRx.QuickView(),Qbeta.View(),P,det); 
	while(N1>0 && QRx.diag()(N1-1)==T(0)) --N1;
      }
#undef AA

      ~QRPDiv() {}

      //
      // Div, DivEq
      //

      template <class T1> void DoLDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(QRx.IsSquare());
	TMVAssert(m.colsize() == QRx.colsize());
	TMVAssert(!istrans);
	QR_LDivEq(QRx,Qbeta,m,N1);
	m /= P; 
      } 

      template <class T1> void DoRDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(QRx.IsSquare());
	TMVAssert(m.rowsize() == QRx.rowsize());
	TMVAssert(!istrans);
	QR_RDivEq(QRx,Qbeta,m%=P,N1);
      }

      template <class T1, class T2> void DoLDiv(const GenMatrix<T1>& m, 
	  const MatrixView<T2>& x) const 
      {
	TMVAssert(m.rowsize() == x.rowsize());
	TMVAssert(m.colsize() == istrans ? QRx.rowsize() : QRx.colsize());
	TMVAssert(x.colsize() == istrans ? QRx.colsize() : QRx.rowsize());
	if (istrans) {
	  QR_RDiv(QRx,Qbeta,m.QuickTranspose()%P,x.QuickTranspose(),N1);
	} else { 
	  QR_LDiv(QRx,Qbeta,m,x,N1);
	  x /= P; 
	} 
      } 

      template <class T1, class T2> void DoRDiv(const GenMatrix<T1>& m, 
	  const MatrixView<T2>& x) const 
      {
	TMVAssert(m.colsize() == x.colsize());
	TMVAssert(m.rowsize() == istrans ? QRx.colsize() : QRx.rowsize());
	TMVAssert(x.rowsize() == istrans ? QRx.rowsize() : QRx.colsize());
	if (istrans) {
	  QR_LDiv(QRx,Qbeta,m.QuickTranspose(),x.QuickTranspose(),N1);
	  x.QuickTranspose() /= P; 
	} else {
	  QR_RDiv(QRx,Qbeta,m%P,x,N1); 
	}
      }

#include "TMV_AuxAllDiv.h"

      //
      // Determinant, Inverse
      //

      T Det() const;
      Matrix<T,ColMajor> Inverse() const;
      Matrix<T,ColMajor> DoInverseATA() const;
      inline bool Singular() const { return Det() == T(0); }
      inline Matrix<T,ColMajor> InverseATA() const
      {
#ifdef TMVDEBUG
	if (istrans) {
	  cout<<"Warning InverseATA called for short matrix in QRPDiv\n";
	  cout<<"The result is really (AAt)^-1\n";
	}
#endif
	return DoInverseATA();
      }

      //
      // Access Decomposition
      //

      Matrix<T,ColMajor> QR_GetQ() const;
      Matrix<T,ColMajor> QR_GetR() const;
      inline Permutation QR_GetP() const { return P; }
      inline bool QR_IsTrans() const { return istrans; }

      inline std::string Type() const
      { return std::string("QRPDiv<") + tmv::Type(T()) + ">"; }

    protected :
      bool istrans;
      Matrix<T,ColMajor> QRx;
      Vector<T> Qbeta;
      Permutation P;
      mutable T det;
      mutable bool donedet;
      size_t N1;
  };


} // namespace mv

#endif
