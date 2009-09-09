//---------------------------------------------------------------------------
//
// This file defines the TMV Divider class.
//
// There are currently 4 algorithms for doing division (and Inverse and Det)
//
// LU Decomposition
// QR Decomposition (with or without Permutation)
// Singular Value Decomposition (compact or full)
// Cholskey (only for SymMatrix)
//
// To tell a Matrix to use a particular algorithm, use the command:
// m.DivideUsing(ALG)
// where ALG is LU, QR, QRP, SV or CH  for the algorithms above.
//
// (There is also SVS which cannot actual perform divisions, but is 
// a quicker SV algorithm that doesn't keep track of the U, V matrices.
// It is useful is for just want the singular values of a matrix.
// Similarly, SVU, SVV only keep S and either U or V.)
//
// The default algorithm is LU for square matrices or QR for non-square.
//
// By default, the appropriate Divider class is created the first
// time it is needed (eg. when the statement v = b/m is called).
// However, you can also setup the Divider class beforehand manually
// by calling m.SetDiv().
//
// You can also query whether the Divider class is already set up.
// This will only be true, if it was previously set up, _and_ the 
// Matrix hasn't been modified since then.
//
// If you want access to the various Divider functions directly,
// They can be accessed by:
//
// m.LUD()
// m.QRD()
// m.SVD()
// m.CHD()
//
// The one of these that is probably most useful to access is SVD(),
// since it is generally a good idea to look for small 
// singular values and zero them out before using SVD for division.
//
// To set to zero all singular value which are less than thresh * 
// the largest singular value use:
//
// m.SVD()->SetThresh(thresh);
//
// To use only the largest nsv singular values use:
//
// m.SVD()->SetTop(nsv);
//
// Also, the singular value decomposition can be used for principal
// component analysis of a Matrix.  The principal component vectors
// are the rows of V.  You can access the decomposition using:
//
// m.SVD()->GetU();
// m.SVD()->GetS(); // The diagonal vector of the (diagonal) matrix W
// m.SVD()->GetV();
//
//


#ifndef TMV_Divider_H
#define TMV_Divider_H

#include "TMV_Matrix.h"

namespace tmv {

#define RT RealType(T)
#define CT ComplexType(T)

  template <class T> class Divider 
  {

    public :

      Divider() {}
      virtual ~Divider() {}

      virtual inline bool IsSV() const { return false; }

      virtual T Det() const =0;
      virtual void Inverse(const MatrixView<T>& minv) const =0;
      virtual void InverseATA(const MatrixView<T>& minv) const =0;
      virtual bool Singular() const =0;
      virtual inline RealType(T) Norm2() const 
      { TMVAssert(FALSE); return RealType(T)(0); }
      virtual inline RealType(T) Condition() const 
      { TMVAssert(FALSE); return RealType(T)(0); }

#define DefDivEq(T) \
      virtual void LDivEq(const MatrixView<T>&) const =0; \
      virtual void RDivEq(const MatrixView<T>&) const =0; \

      DefDivEq(RT)
      DefDivEq(CT)
#undef DefDivEq

#define DefDiv(T1,T2) \
      virtual void LDiv(const GenMatrix<T1>& b, \
	  const MatrixView<T2>& x) const =0; \
      virtual void RDiv(const GenMatrix<T1>& b, \
	  const MatrixView<T2>& x) const =0; \

      DefDiv(RT,RT)
      DefDiv(RT,CT)
      DefDiv(CT,CT)
#undef DefDiv

      virtual string Type() const = 0;
      virtual DivType GetDivType() const = 0;

      virtual bool CheckDecomp(const BaseMatrix<T>& m, ostream* fout) const=0;
  };

#undef RT
#undef CT

  template <class T> inline string Type(const Divider<T>& d)
  { return string("Divider<")+tmv::Type(T())+">{"+d.Type()+"}"; }

} // namespace tmv

#endif
