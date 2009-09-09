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
// quite a bit faster, but it will fail for singular matrices.
// QRPDiv is the version including the column permutations.  
//
// You can control whether QRP does a strict reordering so that the 
// diagonal elements of R are in decreasing order of absolute value 
// (which I call StrictQRP), or whether they are just reordered well
// enough to put the correct zeros at the end (which I call LooseQRP) using
// the global variable tmv::StrictQRP.  The default value is false, which
// is faster, but possibly less accurate for some matrices.
//


#ifndef TMV_QRDiv_H
#define TMV_QRDiv_H

#include "TMV_Vector.h"
#include "TMV_Matrix.h"
#include "TMV_Divider.h"
#include "TMV_TriMatrix.h"

namespace tmv {

  extern bool RecursiveQR;
  extern bool StrictQRP;
  
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
      const MatrixView<T>& QRx, const VectorView<T>& Qbeta, size_t* P, T& det);
  template <class T> void QRP_Decompose(
      const MatrixView<T>& QRx, const VectorView<T>& Qbeta, size_t* P)
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
      const MatrixView<T>& Q, const MatrixView<T>& R, size_t* P, T& det);
  template <class T> void QRP_Decompose(
      const MatrixView<T>& Q, const MatrixView<T>& R, size_t* P)
  { T d=0; QRP_Decompose(Q,R,P,d); }
  // Decompose A (input as Q) into Q R P.
 
  //template <class T> void QR_Decompose(
      //const MatrixView<T>& QRx, vector<GenUpperTriMatrix<T>*>& ZZ, T& det);
  //template <class T> void QR_Decompose(
      //const MatrixView<T>& QRx, vector<GenUpperTriMatrix<T>*>& ZZ)
  //{ T d=0; QR_Decompose(QRx,ZZ,d); }
  // Decompose returning the Z matrices of the BlockHouseholder
  // formulation, rather than just the betas.
  
  template <class T> void GetQFromQR(const MatrixView<T>& Q,
      const GenVector<T>& Qbeta);
  template <class T1, class T2> void Q_LDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& Qbeta,
      const MatrixView<T2>& m);
  template <class T1, class T2> void Q_RDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& Qbeta,
      const MatrixView<T2>& m);
  template <class T1, class T2> inline void Q_LMultEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& Qbeta,
      const MatrixView<T2>& m)
  { Q_RDivEq(Q,Qbeta,m.Adjoint()); }
  template <class T1, class T2> inline void Q_RMultEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& Qbeta,
      const MatrixView<T2>& m)
  { Q_LDivEq(Q,Qbeta,m.Adjoint()); }
  //template <class T1, class T2> void Q_LDivEq(
      //const GenMatrix<T1>& Q, const vector<GenUpperTriMatrix<T1>*> ZZ,
      //const MatrixView<T2>& m);
  //template <class T1, class T2> void Q_RDivEq(
      //const GenMatrix<T1>& Q, const vector<GenUpperTriMatrix<T1>*> ZZ,
      //const MatrixView<T2>& m);

  template <class T> bool QR_DownDate(
      const MatrixView<T>& R, const GenVector<T>& z);
  // Given A = QR, where one of the rows of A is z,
  // find the QR decomposition for A' which is A with the row z removed.
  // This version only updates R, ignoring the Q matrix.

  template <class T> bool QR_DownDate(
      const MatrixView<T>& Q, const MatrixView<T>& R, const GenVector<T>& z);
  // Given A = QR, where the _last_ row of A is z,
  // find the QR decomposition for A' which is A with the row z removed.

  template <class T1, class T2, class T3> void QR_LDiv(
      const GenMatrix<T1>& QR, const GenVector<T1>& Qbeta, const size_t* P,
      const GenMatrix<T2>& m, const MatrixView<T3>& x, size_t N1=0);
  template <class T1, class T2> void QR_LDivEq(
      const GenMatrix<T1>& QR, const GenVector<T1>& Qbeta, const size_t* P,
      const MatrixView<T2>& m, size_t N1=0);
  template <class T1, class T2, class T3> void QR_RDiv(
      const GenMatrix<T1>& QR, const GenVector<T1>& Qbeta, const size_t* P,
      const GenMatrix<T2>& m, const MatrixView<T3>& x, size_t N1=0);
  template <class T1, class T2> void QR_RDivEq(
      const GenMatrix<T1>& QR, const GenVector<T1>& Qbeta, const size_t* P,
      const MatrixView<T2>& m, size_t N1=0);
  //template <class T1, class T2, class T3> void QR_LDiv(
      //const GenMatrix<T1>& QR, const vector<GenUpperTriMatrix<T1>*>& ZZ,
      //const size_t* P, const GenMatrix<T2>& m, const MatrixView<T3>& x,
      //size_t N1=0);
  //template <class T1, class T2> void QR_LDivEq(
      //const GenMatrix<T1>& QR, const vector<GenUpperTriMatrix<T1>*>& ZZ,
      //const size_t* P, const MatrixView<T2>& m, size_t N1=0);
  //template <class T1, class T2, class T3> void QR_RDiv(
      //const GenMatrix<T1>& QR, const vector<GenUpperTriMatrix<T1>*>& ZZ,
      //const size_t* P, const GenMatrix<T2>& m, const MatrixView<T3>& x,
      //size_t N1=0);
  //template <class T1, class T2> void QR_RDivEq(
      //const GenMatrix<T1>& QR, const vector<GenUpperTriMatrix<T1>*>& ZZ,
      //const size_t* P, const MatrixView<T2>& m, size_t N1=0);

#define CT complex<T>
  template <class T> void QR_LDivEq(
      const GenMatrix<CT>& QR, const GenVector<CT>& Qbeta, const size_t* P,
      const MatrixView<T>& m, size_t N1=0)
  { TMVAssert(false); }
  template <class T> void QR_RDivEq(
      const GenMatrix<CT>& QR, const GenVector<CT>& Qbeta, const size_t* P,
      const MatrixView<T>& m, size_t N1=0)
  { TMVAssert(false); }
  //template <class T> void QR_LDivEq(
      //const GenMatrix<CT>& QR, const vector<GenUpperTriMatrix<CT>*>& ZZ,
      //const size_t* P, const MatrixView<T>& m, size_t N1=0)
  //{ TMVAssert(false); }
  //template <class T> void QR_RDivEq(
      //const GenMatrix<CT>& QR, const vector<GenUpperTriMatrix<CT>*>& ZZ, 
      //const size_t* P, const MatrixView<T>& m, size_t N1=0)
  //{ TMVAssert(false); }
#undef CT

  template <class T> class QRDiv : 
    virtual public Divider<T> 
  {

    public :

      //
      // Constructors
      //

      QRDiv(const GenMatrix<T>& A, bool _inplace);
      ~QRDiv();

      //
      // Div, DivEq
      //

      template <class T1> void DoLDivEq(const MatrixView<T1>& m) const
      {
	TMVAssert(QRx.IsSquare());
	TMVAssert(m.colsize() == QRx.colsize());
	if (istrans) QR_LDivEq(QRx,Qbeta,0,m.Transpose(),QRx.rowsize());
	else QR_LDivEq(QRx,Qbeta,0,m,QRx.rowsize());
      } 

      template <class T1> void DoRDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(QRx.IsSquare());
	TMVAssert(m.rowsize() == QRx.rowsize());
	if (istrans) QR_RDivEq(QRx,Qbeta,0,m.Transpose(),QRx.rowsize());
	else QR_RDivEq(QRx,Qbeta,0,m,QRx.rowsize());
      } 

      template <class T1, class T2> void DoLDiv(const GenMatrix<T1>& m, 
	  const MatrixView<T2>& x) const 
      {
	TMVAssert(m.rowsize() == x.rowsize());
	TMVAssert(m.colsize() == istrans ? QRx.rowsize() : QRx.colsize());
	TMVAssert(x.colsize() == istrans ? QRx.colsize() : QRx.rowsize());
	if (istrans) QR_RDiv(QRx,Qbeta,0,m.Transpose(),x.Transpose(),QRx.rowsize());
	else QR_LDiv(QRx,Qbeta,0,m,x,QRx.rowsize());
      } 

      template <class T1, class T2> void DoRDiv(const GenMatrix<T1>& m, 
	  const MatrixView<T2>& x) const 
      {
	TMVAssert(m.colsize() == x.colsize());
	TMVAssert(m.rowsize() == istrans ? QRx.colsize() : QRx.rowsize());
	TMVAssert(x.rowsize() == istrans ? QRx.rowsize() : QRx.colsize());
	if (istrans) QR_LDiv(QRx,Qbeta,0,m.Transpose(),x.Transpose(),QRx.rowsize());
	else QR_RDiv(QRx,Qbeta,0,m,x,QRx.rowsize());
      } 

#include "TMV_AuxAllDiv.h"

      //
      // Determinant, Inverse
      //

      T Det() const;
      void Inverse(const MatrixView<T>& minv) const;
      void DoInverseATA(const MatrixView<T>& minv) const;
      inline bool Singular() const { return Det() == T(0); }
      inline void InverseATA(const MatrixView<T>& minv) const
      {
#ifdef TMVDEBUG
	if (istrans) {
	  cout<<"Warning InverseATA called for short matrix in QRDiv\n";
	  cout<<"The result is really (AAt)^-1\n";
	}
#endif
	DoInverseATA(minv);
      }


      //
      // Access Decomposition
      //

      inline bool IsTrans() const { return istrans; }
      Matrix<T> GetQ() const;
      ConstUpperTriMatrixView<T> GetR() const
      { return UpperTriMatrixViewOf(QRx,NonUnitDiag); }
      const GenMatrix<T>& GetQR() const { return QRx; }
      const GenVector<T>& GetQbeta() const { return Qbeta; }

      inline std::string Type() const
      { return std::string("QRDiv<") + tmv::Type(T()) + ">"; }
      inline DivType GetDivType() const { return QR; }

      bool CheckDecomp(const BaseMatrix<T>& m, ostream* fout) const;

    protected :
      bool istrans;
      bool inplace;
      T* Aptr;
      MatrixView<T> QRx;
      Vector<T> Qbeta;
      mutable T det;
      mutable bool donedet;
  };

  template <class T> class QRPDiv : 
    virtual public Divider<T> 
  {

    public :

      //
      // Constructors
      //

      QRPDiv(const GenMatrix<T>& A, bool _inplace);
      ~QRPDiv();

      //
      // Div, DivEq
      //

      template <class T1> void DoLDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(QRx.IsSquare());
	TMVAssert(m.colsize() == QRx.colsize());
	if (istrans) QR_LDivEq(QRx,Qbeta,P,m.Transpose(),N1);
	else QR_LDivEq(QRx,Qbeta,P,m,N1);
      } 

      template <class T1> void DoRDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(QRx.IsSquare());
	TMVAssert(m.rowsize() == QRx.rowsize());
	if (istrans) QR_RDivEq(QRx,Qbeta,P,m.Transpose(),N1);
	else QR_RDivEq(QRx,Qbeta,P,m,N1);
      }

      template <class T1, class T2> void DoLDiv(const GenMatrix<T1>& m, 
	  const MatrixView<T2>& x) const 
      {
	TMVAssert(m.rowsize() == x.rowsize());
	TMVAssert(m.colsize() == istrans ? QRx.rowsize() : QRx.colsize());
	TMVAssert(x.colsize() == istrans ? QRx.colsize() : QRx.rowsize());
	if (istrans) QR_RDiv(QRx,Qbeta,P,m.Transpose(),x.Transpose(),N1);
	else QR_LDiv(QRx,Qbeta,P,m,x,N1);
      } 

      template <class T1, class T2> void DoRDiv(const GenMatrix<T1>& m, 
	  const MatrixView<T2>& x) const 
      {
	TMVAssert(m.colsize() == x.colsize());
	TMVAssert(m.rowsize() == istrans ? QRx.colsize() : QRx.rowsize());
	TMVAssert(x.rowsize() == istrans ? QRx.rowsize() : QRx.colsize());
	if (istrans) QR_LDiv(QRx,Qbeta,P,m.Transpose(),x.Transpose(),N1);
	else QR_RDiv(QRx,Qbeta,P,m,x,N1); 
      }

#include "TMV_AuxAllDiv.h"

      //
      // Determinant, Inverse
      //

      T Det() const;
      void Inverse(const MatrixView<T>& minv) const;
      void DoInverseATA(const MatrixView<T>& minv) const;
      inline bool Singular() const { return Det() == T(0); }
      inline void InverseATA(const MatrixView<T>& minv) const
      {
#ifdef TMVDEBUG
	if (istrans) {
	  cout<<"Warning InverseATA called for short matrix in QRPDiv\n";
	  cout<<"The result is really (AAt)^-1\n";
	}
#endif
	DoInverseATA(minv);
      }

      //
      // Access Decomposition
      //

      inline bool IsTrans() const { return istrans; }
      Matrix<T> GetQ() const;
      inline ConstUpperTriMatrixView<T> GetR() const
      { return UpperTriMatrixViewOf(QRx); }
      inline const GenMatrix<T>& GetQRx() const { return QRx; }
      inline const GenVector<T>& GetQbeta() const { return Qbeta; }
      inline const size_t* GetP() const { return P; }

      inline std::string Type() const
      { return std::string("QRPDiv<") + tmv::Type(T()) + ">"; }
      inline DivType GetDivType() const { return QRP; }

      bool CheckDecomp(const BaseMatrix<T>& m, ostream* fout) const;

    protected :
      bool istrans;
      bool inplace;
      T* Aptr;
      MatrixView<T> QRx;
      Vector<T> Qbeta;
      size_t* P;
      mutable T det;
      mutable bool donedet;
      size_t N1;
  };

} // namespace mv

#endif
