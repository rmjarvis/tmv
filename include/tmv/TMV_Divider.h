///////////////////////////////////////////////////////////////////////////////
// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2009                                                 //
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

#include "tmv/TMV_BaseMatrix.h"

namespace tmv {

#define RT RealType(T)
#define CT ComplexType(T)

  template <class T> 
  class Divider 
  {

  public :

    Divider() {}
    virtual ~Divider() {}

    virtual inline bool IsSV() const { return false; }

    virtual T Det() const =0;
    virtual RealType(T) LogDet(T* sign) const =0;
    virtual void InverseATA(const MatrixView<T>& minv) const =0;
    virtual bool Singular() const =0;
    virtual inline RealType(T) Norm2() const 
    { TMVAssert(FALSE); return RealType(T)(0); }
    virtual inline RealType(T) Condition() const 
    { TMVAssert(FALSE); return RealType(T)(0); }

#define DefDivEq(T) \
    virtual void LDivEq(const MatrixView<T>&) const =0; \
    virtual void RDivEq(const MatrixView<T>&) const =0; \
    virtual void Inverse(const MatrixView<T>& minv) const =0; \

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

    virtual bool CheckDecomp(const BaseMatrix<T>& m, 
        std::ostream* fout) const=0;
  };

#undef RT
#undef CT

  template <class T> 
  inline std::string TypeText(const Divider<T>& d)
  { return std::string("Divider<")+TypeText(T())+">"; }

} // namespace tmv

#endif
