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
#ifndef TMV_SVD_H
#define TMV_SVD_H

#include "tmv/TMV_Divider.h"
#include "tmv/TMV_BaseDiagMatrix.h"

namespace tmv {

  // Decompose A (input as U) into A = U S V
  // where U is column-unitary, S is diagonal, and V is unitary.
  // U must have U.nrows() >= U.ncols().
  // S must have S.size() == U.ncols().
  // V must have V.nrows() == V.ncols() == U.ncols().
  // If StoreU is false, then U will be junk on output.
  template <class T> 
  void SV_Decompose(const MatrixView<T>& U,
      const DiagMatrixView<RealType(T)>& S, 
      const MatrixView<T>& V, bool StoreU=true);

  // The same, but don't return V:
  template <class T> 
  void SV_Decompose(const MatrixView<T>& U,
      const DiagMatrixView<RealType(T)>& S, bool StoreU);

  template <class T> 
  class SVDiv : 
    public Divider<T> 
  {

  public :

    //
    // Constructors
    //

    SVDiv(const GenMatrix<T>& A, bool _inplace);
    ~SVDiv();

    //
    // Div, DivEq
    //

    template <class T1> 
    void DoLDivEq(const MatrixView<T1>& m) const;
    template <class T1> 
    void DoRDivEq(const MatrixView<T1>& m) const;
    template <class T1, class T2> 
    void DoLDiv(const GenMatrix<T1>& m, const MatrixView<T2>& x) const;
    template <class T1, class T2> 
    void DoRDiv(const GenMatrix<T1>& m, const MatrixView<T2>& x) const;

#include "tmv/TMV_AuxAllDiv.h"

    //
    // Determinant, Inverse
    //

    T Det() const;
    RealType(T) LogDet(T* sign) const;
    template <class T1> 
    void DoInverse(const MatrixView<T1>& minv) const;
    void DoInverseATA(const MatrixView<T>& minv) const;
    bool Singular() const;
    RealType(T) Norm2() const;
    RealType(T) Condition() const;

    //
    // Determine which (if any) S values to zero out
    //

    void Thresh(RealType(T) toler, std::ostream* debugout=0) const;
    void Top(int neigen, std::ostream* debugout=0) const;
    int GetKMax() const;

    //
    // Access Decomposition
    //

    ConstMatrixView<T> GetU() const;
    ConstDiagMatrixView<RealType(T)> GetS() const;
    ConstMatrixView<T> GetV() const;

    bool CheckDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

  protected:

    struct SVDiv_Impl;
    std::auto_ptr<SVDiv_Impl> pimpl;

    size_t colsize() const;
    size_t rowsize() const;

  private :

    inline SVDiv(const SVDiv<T>&) : pimpl(0)
    { TMVAssert(FALSE); }
    inline SVDiv<T>& operator=(const SVDiv<T>&)
    { TMVAssert(FALSE); return *this; }

  };

} // namespace tmv;

#endif
