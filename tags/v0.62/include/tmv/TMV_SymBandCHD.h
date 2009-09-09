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


#ifndef TMV_SymBandCHD_H
#define TMV_SymBandCHD_H

#include "tmv/TMV_SymDivider.h"
#include "tmv/TMV_BaseSymBandMatrix.h"
#include "tmv/TMV_BaseDiagMatrix.h"

namespace tmv {

  // Decompose A into L*Lt.
  // On output, L = A.LowerBand().
  template <class T> 
  void CH_Decompose(const SymBandMatrixView<T>& A);

  // Decompose Tridiagonal A into L D Lt
  template <class T> 
  void LDL_Decompose(const SymBandMatrixView<T>& A);

  template <class T> 
  class HermBandCHDiv : 
    public SymDivider<T> 
  {

  public :

    //
    // Constructors
    //

    HermBandCHDiv(const GenSymBandMatrix<T>& A, bool _inplace);
    ~HermBandCHDiv();

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

    const BandMatrix<T> GetL() const;
    const DiagMatrix<T> GetD() const;
    const GenSymBandMatrix<T>& GetLL() const;

    bool CheckDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

  private :
    struct HermBandCHDiv_Impl;
    std::auto_ptr<HermBandCHDiv_Impl> pimpl;

    size_t colsize() const;
    size_t rowsize() const;

  private :

    inline HermBandCHDiv(const HermBandCHDiv<T>&) : pimpl(0)
    { TMVAssert(FALSE); }
    inline HermBandCHDiv<T>& operator=(const HermBandCHDiv<T>&)
    { TMVAssert(FALSE); return *this; }

  };

} // namespace tmv

#endif
