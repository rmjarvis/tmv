//---------------------------------------------------------------------------
//
// This file contains the code for doing division of BandMatrices using 
// QR Decomposition.
//
// I don't implement the QRP method, since the permutations screw up the
// band structure.  It could be implemented using a similar technique 
// as I used for the BandLUDiv class, but if you want to use it you need 
// to cast the band matrix to a regular matrix first.
//


#ifndef TMV_BandQRDiv_H
#define TMV_BandQRDiv_H

#include "TMV_BandMatrix.h"
#include "TMV_Divider.h"
#include "TMV_QRDiv.h"

namespace tmv {

  template <class T> void BandQR_Decompose(
      const BandMatrixView<T>& QRx, const VectorView<T>& Qbeta, T& det);
  template <class T> inline void BandQR_Decompose(
      const BandMatrixView<T>& QRx, const VectorView<T>& Qbeta)
  { T d; BandQR_Decompose(QRx,Qbeta,d); }
  // Decompose A (input as QRx) into Q R.
  // On output, Q is stored in the lower triangle part of QRx as
  // Householder vectors.  The upper triangle part and Qbeta hold R.

  template <class T1, class T2> void BandQR_LDivEq(
      const GenBandMatrix<T1>& QRx, const GenVector<T1>& Qbeta, 
      const MatrixView<T2>& m);
  template <class T1, class T2> void BandQR_RDivEq(
      const GenBandMatrix<T1>& QRx, const GenVector<T1>& Qbeta, 
      const MatrixView<T2>& m);

  template <class T1, class T2, class T3> void BandQR_LDiv(
	const GenBandMatrix<T1>& QR, const GenVector<T1>& Qbeta,
	const GenMatrix<T2>& m, const MatrixView<T3>& x);
  template <class T1, class T2, class T3> void BandQR_RDiv(
	const GenBandMatrix<T1>& QR, const GenVector<T1>& Qbeta, 
	const GenMatrix<T2>& m, const MatrixView<T3>& x);

#define CT complex<T>
  template <class T> inline void BandQR_LDivEq(
      const GenBandMatrix<CT>& , const GenVector<CT>& ,
      const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void BandQR_RDivEq(
      const GenBandMatrix<CT>& , const GenVector<CT>& ,
      const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void BandQR_LDiv(
      const GenBandMatrix<CT>& , const GenVector<CT>& ,
      const GenMatrix<T>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void BandQR_RDiv(
      const GenBandMatrix<CT>& , const GenVector<CT>& ,
      const GenMatrix<T>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
#undef CT

  template <class T1, class T2> void BandQ_LDivEq(
      const GenBandMatrix<T1>& Q, const GenVector<T1>& Qbeta, 
      const MatrixView<T2>& m);
  template <class T1, class T2> void BandQ_RDivEq(
      const GenBandMatrix<T1>& Q, const GenVector<T1>& Qbeta, 
      const MatrixView<T2>& m);
  template <class T> Matrix<T> GetQFromBandQR(
      const GenBandMatrix<T>& QRx, const GenVector<T>& Qbeta);

  template <class T> class BandQRDiv : 
    public Divider<T> 
  {

    public :

      //
      // Constructors
      //

      BandQRDiv(const GenBandMatrix<T>& A, bool _inplace);
      ~BandQRDiv() {}

      //
      // Div, DivEq
      //
      
      template <class T1> inline void DoLDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(QRx.IsSquare());
	TMVAssert(m.colsize() == QRx.colsize());
	if (istrans) BandQR_RDivEq(QRx,Qbeta,m.Transpose());
	else BandQR_LDivEq(QRx,Qbeta,m);
      }

      template <class T1> inline void DoRDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(QRx.IsSquare());
	TMVAssert(m.rowsize() == QRx.colsize());
	if (istrans) BandQR_LDivEq(QRx,Qbeta,m.Transpose());
	else BandQR_RDivEq(QRx,Qbeta,m);
      }

      template <class T1, class T2> inline void DoLDiv(
	  const GenMatrix<T1>& m, const MatrixView<T2>& x) const 
      { 
	TMVAssert(m.rowsize() == x.rowsize());
	TMVAssert(m.colsize() == istrans ? QRx.rowsize() : QRx.colsize());
	TMVAssert(x.colsize() == istrans ? QRx.colsize() : QRx.rowsize());
	if (istrans) BandQR_RDiv(QRx,Qbeta,m.Transpose(),x.Transpose()); 
	else BandQR_LDiv(QRx,Qbeta,m,x); 
      }

      template <class T1, class T2> inline void DoRDiv(
	  const GenMatrix<T1>& m, const MatrixView<T2>& x) const 
      { 
	TMVAssert(m.colsize() == x.colsize());
	TMVAssert(m.rowsize() == istrans ? QRx.colsize() : QRx.rowsize());
	TMVAssert(x.rowsize() == istrans ? QRx.rowsize() : QRx.colsize());
	if (istrans) BandQR_LDiv(QRx,Qbeta,m.Transpose(),x.Transpose()); 
	else BandQR_RDiv(QRx,Qbeta,m,x); 
      }

#include "TMV_AuxAllDiv.h"

      //
      // Determinant, Inverse
      //

      T Det() const;
      void Inverse(const MatrixView<T>& minv) const;
      void DoInverseATA(const MatrixView<T>& minv) const;
      inline void InverseATA(const MatrixView<T>& minv) const
      {
#ifdef TMVDEBUG
	if (istrans) {
	  cout<<"Warning InverseATA called for short matrix in BandQRDiv\n";
	  cout<<"The result is really (AAt)^-1\n";
	}
#endif
	DoInverseATA(minv);
      }
      inline bool Singular() const { return Det() == T(0); }

      //
      // Access Decomposition
      //

      inline bool IsTrans() const { return istrans; }
      Matrix<T> GetQ() const;
      inline ConstBandMatrixView<T> GetR() const 
      { return QRx.SubBandMatrix(0,QRx.rowsize(),0,QRx.rowsize(),0,QRx.nhi()); }
      inline const GenBandMatrix<T>& GetQR() const { return QRx; }
      inline const GenVector<T>& GetQbeta() const { return Qbeta; }

      inline std::string Type() const
      { return std::string("BandQRDiv<") + tmv::Type(T()) + ">"; }
      inline DivType GetDivType() const { return QR; }

      bool CheckDecomp(const BaseMatrix<T>& m, ostream* fout) const;

    private :

      bool istrans;
      bool inplace;
      auto_array<T> Aptr1;
      T* Aptr;
      BandMatrixView<T> QRx;
      Vector<T> Qbeta;
      mutable T det;
      mutable bool donedet;
  };

} // namespace mv

#endif
