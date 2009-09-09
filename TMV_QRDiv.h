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
// If R is singular, then you need QRP Decomposition (see TMV_QRPDiv.h).
//


#ifndef TMV_QRDiv_H
#define TMV_QRDiv_H

#include "TMV_Vector.h"
#include "TMV_Matrix.h"
#include "TMV_Divider.h"
#include "TMV_TriMatrix.h"

namespace tmv {

  template <class T> void QR_Decompose(
      const MatrixView<T>& QRx, const VectorView<T>& beta, T& det);
  template <class T> inline void QR_Decompose(
      const MatrixView<T>& QRx, const VectorView<T>& beta)
  { T d=0; QR_Decompose(QRx,beta,d); }
  // Decompose A (input as QRx) into Q R.
  // On output, Q is stored in the lower triangle part of QRx as
  // Householder vectors, and the beta vector.
  // R is the upper triangle part of QRx
 
  template <class T> void QR_Decompose(
      const MatrixView<T>& Q, const UpperTriMatrixView<T>& R, T& det);
  template <class T> inline void QR_Decompose(
      const MatrixView<T>& Q, const UpperTriMatrixView<T>& R)
  { T d=0; QR_Decompose(Q,R,d); }
  // Decompose A (input as Q) into Q R.
 
  template <class T> void GetQFromQR(const MatrixView<T>& Q,
      const GenVector<T>& beta);
  template <class T, class T1> void Q_LDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& beta,
      const MatrixView<T>& m);
  template <class T, class T1> void Q_RDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& beta,
      const MatrixView<T>& m);
  template <class T, class T1> inline void Q_LMultEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& beta,
      const MatrixView<T>& m)
  { Q_RDivEq(Q,beta,m.Adjoint()); }
  template <class T, class T1> inline void Q_RMultEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& beta,
      const MatrixView<T>& m)
  { Q_LDivEq(Q,beta,m.Adjoint()); }

  template <class T> void QR_Update(
      const UpperTriMatrixView<T>& R, const MatrixView<T>& A);
  // Given that A0 = Q0 R0
  // Find R1, so that [ A0 ] = Q1 R1
  //                  [ A  ] 
  // On input R is R0; on output it is R1.
  template <class T> inline void QR_Update(
      const UpperTriMatrixView<T>& R, const VectorView<T>& z)
  { return QR_Update(R,RowVectorViewOf(z)); }
  template <class T> bool QR_Downdate(
      const UpperTriMatrixView<T>& R, const MatrixView<T>& A);
  // The opposite of an Update:
  // Given that [ A0 ] = Q1 R1
  //            [ A  ]
  // Find R0, so that A0 = Q0 R0
  // On input R is R1; on output it is R0.
  template <class T> inline bool QR_Downdate(
      const UpperTriMatrixView<T>& R, const VectorView<T>& z)
  { return QR_Downdate(R,RowVectorViewOf(z)); }

  template <class T, class T1, class T2> void QR_LDiv(
      const GenMatrix<T1>& QR, const GenVector<T1>& beta, const size_t* P,
      const GenMatrix<T2>& m, const MatrixView<T>& x, size_t N1=0);
  template <class T, class T1> void QR_LDivEq(
      const GenMatrix<T1>& QR, const GenVector<T1>& beta, const size_t* P,
      const MatrixView<T>& m, size_t N1=0);
  template <class T, class T1, class T2> void QR_RDiv(
      const GenMatrix<T1>& QR, const GenVector<T1>& beta, const size_t* P,
      const GenMatrix<T2>& m, const MatrixView<T>& x, size_t N1=0);
  template <class T, class T1> void QR_RDivEq(
      const GenMatrix<T1>& QR, const GenVector<T1>& beta, const size_t* P,
      const MatrixView<T>& m, size_t N1=0);

#define CT complex<T>
  template <class T> inline void QR_LDivEq(
      const GenMatrix<CT>& , const GenVector<CT>& , const size_t* ,
      const MatrixView<T>& , size_t =0)
  { TMVAssert(FALSE); }
  template <class T> inline void QR_RDivEq(
      const GenMatrix<CT>& , const GenVector<CT>& , const size_t* ,
      const MatrixView<T>& , size_t =0)
  { TMVAssert(FALSE); }
  template <class T> inline void QR_LDiv(
      const GenMatrix<CT>& , const GenVector<CT>& , const size_t* ,
      const GenMatrix<T>& , const MatrixView<T>& , size_t =0)
  { TMVAssert(FALSE); }
  template <class T> inline void QR_RDiv(
      const GenMatrix<CT>& , const GenVector<CT>& , const size_t* ,
      const GenMatrix<T>& , const MatrixView<T>& , size_t =0)
  { TMVAssert(FALSE); }
#undef CT

  template <class T> class QRDiv : 
    public Divider<T> 
  {

    public :

      //
      // Constructors
      //

      QRDiv(const GenMatrix<T>& A, bool _inplace);
      ~QRDiv() {}

      //
      // Div, DivEq
      //

      template <class T1> void DoLDivEq(const MatrixView<T1>& m) const
      {
	TMVAssert(QRx.IsSquare());
	TMVAssert(m.colsize() == QRx.colsize());
	if (istrans) QR_LDivEq(QRx,beta,0,m.Transpose(),QRx.rowsize());
	else QR_LDivEq(QRx,beta,0,m,QRx.rowsize());
      } 

      template <class T1> void DoRDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(QRx.IsSquare());
	TMVAssert(m.rowsize() == QRx.rowsize());
	if (istrans) QR_RDivEq(QRx,beta,0,m.Transpose(),QRx.rowsize());
	else QR_RDivEq(QRx,beta,0,m,QRx.rowsize());
      } 

      template <class T1, class T2> void DoLDiv(const GenMatrix<T1>& m, 
	  const MatrixView<T2>& x) const 
      {
	TMVAssert(m.rowsize() == x.rowsize());
	TMVAssert(m.colsize() == istrans ? QRx.rowsize() : QRx.colsize());
	TMVAssert(x.colsize() == istrans ? QRx.colsize() : QRx.rowsize());
	if (istrans) QR_RDiv(QRx,beta,0,m.Transpose(),x.Transpose(),QRx.rowsize());
	else QR_LDiv(QRx,beta,0,m,x,QRx.rowsize());
      } 

      template <class T1, class T2> void DoRDiv(const GenMatrix<T1>& m, 
	  const MatrixView<T2>& x) const 
      {
	TMVAssert(m.colsize() == x.colsize());
	TMVAssert(m.rowsize() == istrans ? QRx.colsize() : QRx.rowsize());
	TMVAssert(x.rowsize() == istrans ? QRx.rowsize() : QRx.colsize());
	if (istrans) QR_LDiv(QRx,beta,0,m.Transpose(),x.Transpose(),QRx.rowsize());
	else QR_RDiv(QRx,beta,0,m,x,QRx.rowsize());
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
      const GenVector<T>& Getbeta() const { return beta; }

      inline string Type() const
      { return string("QRDiv<") + tmv::Type(T()) + ">"; }
      inline DivType GetDivType() const { return QR; }

      bool CheckDecomp(const BaseMatrix<T>& m, ostream* fout) const;

    protected :
      bool istrans;
      bool inplace;
      auto_array<T> Aptr1;
      T* Aptr;
      MatrixView<T> QRx;
      Vector<T> beta;
      mutable T det;
      mutable bool donedet;
  };

} // namespace mv

#endif
