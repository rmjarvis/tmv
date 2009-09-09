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


//---------------------------------------------------------------------------
//
// This file contains the code for doing division using 
// QRP Decomposition.
//
// In a QR Decomposition, if R is singular, then back-substitution will 
// fail in solving R x = Qt b.  However, there is a trick that still 
// gives us a valid least squares solution to A x = b.  
// This is to use column pivoting to obtain a decomposition:
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
// You can control whether QRP does a strict reordering so that the 
// diagonal elements of R are in decreasing order of absolute value 
// (which I call StrictQRP), or whether they are just reordered well
// enough to put the correct zeros at the end (which I call LooseQRP) using
// the global variable tmv::StrictQRP.  The default value is false, which
// is faster, but possibly less accurate for some matrices.
//


#ifndef TMV_QRPD_H
#define TMV_QRPD_H

#include "TMV_Divider.h"
#include "TMV_BaseTriMatrix.h"

namespace tmv {

  // Decompose A (input as Q) into Q R P.
  template <class T> void QRP_Decompose(
      const MatrixView<T>& Q, const UpperTriMatrixView<T>& R, int* P,
      bool strict=false);
  // Decompose A into Q R P, but don't return Q or P.
  // R is returned as UpperTriMatrixViewOf(A).
  template <class T> void QRP_Decompose(const MatrixView<T>& A,
      bool strict=false);
 
  template <class T> class QRPDiv : 
    public Divider<T> 
  {

    public :

      static bool StrictQRP;

      //
      // Constructors
      //

      QRPDiv(const GenMatrix<T>& A, bool _inplace);
      ~QRPDiv();

      //
      // Div, DivEq
      //

      template <class T1> void DoLDivEq(const MatrixView<T1>& m) const;
      template <class T1> void DoRDivEq(const MatrixView<T1>& m) const;
      template <class T1, class T2> void DoLDiv(const GenMatrix<T1>& m, 
	  const MatrixView<T2>& x) const;
      template <class T1, class T2> void DoRDiv(const GenMatrix<T1>& m, 
	  const MatrixView<T2>& x) const;

      //
      // Determinant, Inverse
      //

      T Det() const;
      RealType(T) LogDet(T* sign) const;
      template <class T1> void DoInverse(const MatrixView<T1>& minv) const;
      void DoInverseATA(const MatrixView<T>& minv) const;
      bool Singular() const;

#include "TMV_AuxAllDiv.h"

      //
      // Access Decomposition
      //

      bool IsTrans() const;
      Matrix<T> GetQ() const;
      ConstUpperTriMatrixView<T> GetR() const;
      const GenMatrix<T>& GetQRx() const;
      const GenVector<T>& Getbeta() const;
      const int* GetP() const;

      bool CheckDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

    protected :

      struct QRPDiv_Impl;
      QRPDiv_Impl* pimpl;

      size_t colsize() const;
      size_t rowsize() const;

    private :

      inline QRPDiv(const QRPDiv<T>&) : pimpl(0)
      { TMVAssert(FALSE); }
      inline QRPDiv<T>& operator=(const QRPDiv<T>&)
      { TMVAssert(FALSE); return *this; }

  };

} // namespace mv

#endif
