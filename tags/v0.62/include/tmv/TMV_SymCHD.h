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
// This file contains the code for doing division using 
// Cholesky Decomposition.
// 
// The algorithm is much like the LU decomposition, but we don't do
// any pivoting, and since the source matrix is symmetric, L = LT
// (or for Hermition, L = Lt).
//


#ifndef TMV_SymCHD_H
#define TMV_SymCHD_H

#include "tmv/TMV_SymDivider.h"
#include "tmv/TMV_BaseTriMatrix.h"
#include "tmv/TMV_BaseBandMatrix.h"

namespace tmv {

  // Decompose A into L*Lt
  // On output L = A.LowerTri()
  template <class T> 
  void CH_Decompose(const SymMatrixView<T>& A);

  // Decompose A into U*P
  // where U is unitary and P is positive definite
  // On ouput U = A
  template <class T> 
  void Polar_Decompose(const MatrixView<T>& A, const SymMatrixView<T>& P);

  template <class T> 
  void Polar_Decompose(const GenBandMatrix<T>& A,
      const MatrixView<T>& U, const SymMatrixView<T>& P);

  template <class T> 
  void SquareRoot(const SymMatrixView<T>& A);

  template <class T> 
  class HermCHDiv : 
    public SymDivider<T> 
  {

  public :

    //
    // Constructors
    //

    HermCHDiv(const GenSymMatrix<T>& A, bool _inplace);
    ~HermCHDiv();

    //
    // Divider Versions of DivEq and Div
    //

    template <class T1> 
    void DoLDivEq(const MatrixView<T1>& m) const;

    template <class T1> 
    void DoRDivEq(const MatrixView<T1>& m) const;

    template <class T1, class T2> 
    void DoLDiv(const GenMatrix<T1>& m1, const MatrixView<T2>& m0) const;

    template <class T1, class T2> 
    void DoRDiv(const GenMatrix<T1>& m1, const MatrixView<T2>& m0) const;

    //
    // Determinant, Inverse
    //

    T Det() const;
    RealType(T) LogDet(T* sign) const;
    template <class T1> 
    void DoInverse(const MatrixView<T1>& minv) const;

    template <class T1> 
    void DoInverse(const SymMatrixView<T1>& minv) const;

    inline void Inverse(const SymMatrixView<RealType(T)>& sinv) const
    {
      TMVAssert(IsReal(T()));
      TMVAssert(sinv.size() == colsize());
      DoInverse(sinv);
    }
    inline void Inverse(const SymMatrixView<ComplexType(T)>& sinv) const
    {
      TMVAssert(sinv.size() == colsize());
      TMVAssert(sinv.isherm());
      DoInverse(sinv);
    }
    void DoInverseATA(const MatrixView<T>& minv) const;
    bool Singular() const;

#include "tmv/TMV_AuxAllDiv.h"

    //
    // Access Decomposition
    //

    const ConstLowerTriMatrixView<T> GetL() const;
    const GenSymMatrix<T>& GetLL() const;

    bool CheckDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

  private :
    struct HermCHDiv_Impl;
    std::auto_ptr<HermCHDiv_Impl> pimpl;

    size_t colsize() const;
    size_t rowsize() const;

  private :

    inline HermCHDiv(const HermCHDiv<T>&) : pimpl(0)
    { TMVAssert(FALSE); }
    inline HermCHDiv<T>& operator=(const HermCHDiv<T>&)
    { TMVAssert(FALSE); return *this; }

  }; // HermCHDiv

} // namespace tmv

#endif
