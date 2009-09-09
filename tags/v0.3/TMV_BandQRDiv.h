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
  template <class T> void BandQR_Decompose(
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
      const GenBandMatrix<CT>& QRx, const GenVector<CT>& Qbeta,
      const MatrixView<T>& m)
  { TMVAssert(false); }
  template <class T> inline void BandQR_RDivEq(
      const GenBandMatrix<CT>& QRx, const GenVector<CT>& Qbeta,
      const MatrixView<T>& m)
  { TMVAssert(false); }
#undef CT

  template <class T> class BandQRDiv : 
    virtual public BaseQRDiv<T> 
  {

    public :

      //
      // Constructors
      //

      explicit BandQRDiv(const GenBandMatrix<T>& A) :
	istrans(A.colsize()<A.rowsize() || (A.IsSquare() && A.nhi()<A.nlo())),
	QRx(istrans?A.rowsize():A.colsize(),
	    istrans?A.colsize():A.rowsize(),
	    istrans?A.nhi():A.nlo(), 
	    min(A.nlo()+A.nhi(),int(istrans?A.colsize():A.rowsize())-1),
	    istrans?TransOf(BaseStorOf(A)):BaseStorOf(A)),
	Qbeta(QRx.rowsize()), det(T(1)), donedet(false)
      {
	if (istrans) {
	  BandMatrixViewOf(QRx.QuickTranspose(),A.nlo(),A.nhi()) = A;
	  for(int i=A.nlo()+1;i<=QRx.nhi();++i) QRx.diag(i).Zero();
	} else {
	  BandMatrixViewOf(QRx,A.nlo(),A.nhi()) = A;
	  for(int i=A.nhi()+1;i<=QRx.nhi();++i) QRx.diag(i).Zero();
	}
	BandQR_Decompose(QRx.QuickView(),Qbeta.View(),det);
      }

      ~BandQRDiv() {}

      //
      // Div, DivEq
      //
      
      template <class T1> inline void DoLDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(QRx.IsSquare());
	if (istrans) {
	  TMVAssert(m.colsize() == QRx.colsize());
	  BandQR_RDivEq(QRx,Qbeta,m.QuickTranspose());
	} else {
	  TMVAssert(m.colsize() == QRx.colsize());
	  BandQR_LDivEq(QRx,Qbeta,m);
	}
      }

      template <class T1> inline void DoRDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(QRx.IsSquare());
	if (istrans) {
	  TMVAssert(m.rowsize() == QRx.colsize());
	  BandQR_LDivEq(QRx,Qbeta,m.QuickTranspose());
	} else {
	  TMVAssert(m.rowsize() == QRx.colsize());
	  BandQR_RDivEq(QRx,Qbeta,m);
	}
      }

      template <class T1, class T2> inline void DoLDiv(
	    const GenMatrix<T1>& m, const MatrixView<T2>& x) const 
	{ 
	  TMVAssert(m.rowsize() == x.rowsize());
	  if (istrans) {
	    TMVAssert(m.colsize() == QRx.rowsize());
	    TMVAssert(x.colsize() == QRx.colsize());
	    BandQR_RDiv(QRx,Qbeta,m.QuickTranspose(),x.QuickTranspose()); 
	  } else {
	    TMVAssert(m.colsize() == QRx.colsize());
	    TMVAssert(x.colsize() == QRx.rowsize());
	    BandQR_LDiv(QRx,Qbeta,m,x); 
	  }
	}

      template <class T1, class T2> inline void DoRDiv(
	    const GenMatrix<T1>& m, const MatrixView<T2>& x) const 
	{ 
	  TMVAssert(m.colsize() == x.colsize());
	  if (istrans) {
	    TMVAssert(m.rowsize() == QRx.colsize());
	    TMVAssert(x.rowsize() == QRx.rowsize());
	    BandQR_LDiv(QRx,Qbeta,m.QuickTranspose(),x.QuickTranspose()); 
	  } else {
	    TMVAssert(m.rowsize() == QRx.rowsize());
	    TMVAssert(x.rowsize() == QRx.colsize());
	    BandQR_RDiv(QRx,Qbeta,m,x); 
	  }
	}

#include "TMV_AuxAllDiv.h"

      //
      // Determinant, Inverse
      //

      T Det() const;
      Matrix<T,ColMajor> Inverse() const;
      Matrix<T,ColMajor> DoInverseATA() const;
      inline Matrix<T,ColMajor> InverseATA() const
      {
#ifdef TMVDEBUG
	if (istrans) {
	  cout<<"Warning InverseATA called for short matrix in BandQRDiv\n";
	  cout<<"The result is really (AAt)^-1\n";
	}
#endif
	return DoInverseATA();
      }
      inline bool Singular() const { return Det() == T(0); }

      //
      // Access Decomposition
      //


      Matrix<T,ColMajor> QR_GetQ() const;
      inline Matrix<T,ColMajor> QR_GetR() const 
      { return Matrix<T,ColMajor>(GetBandR()); }
      inline Permutation QR_GetP() const { return Permutation(Qbeta.size()); }
      inline bool QR_IsTrans() const { return istrans; }

      ConstBandMatrixView<T> GetBandR() const 
      { return QRx.SubBandMatrix(0,QRx.rowsize(),0,QRx.rowsize(),0,QRx.nhi()); }
      const BandMatrix<T,ColMajor>& GetBandQR() const { return QRx; }
      const Vector<T>& GetQbeta() const { return Qbeta; }

      inline std::string Type() const
      { return std::string("BandQRDiv<") + tmv::Type(T()) + ">"; }

    private :

      bool istrans;
      BandMatrix<T,ColMajor> QRx;
      Vector<T> Qbeta;
      mutable T det;
      mutable bool donedet;
  };

} // namespace mv

#endif
