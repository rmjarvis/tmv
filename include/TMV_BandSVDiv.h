///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2007                                                        //
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
// along with this program in the file gpl.txt.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


//---------------------------------------------------------------------------
#ifndef TMV_BandSVDiv_H
#define TMV_BandSVDiv_H

#include "TMV_Divider.h"
//#include "TMV_SVDiv.h"

namespace tmv {

  template <class T> void BandSV_Decompose(
      const GenBandMatrix<T>& A,
      const MatrixView<T>* U, const VectorView<RealType(T)>& S, 
      const MatrixView<T>* V, T& det);
  // Decompose A into U S V
  // where S is a diagonal real matrix, and U,V are unitary matrices.
  // U,S,V are N x N

  template <class T> class BandSVDiv : 
    public Divider<T> 
  {

    public :

      //
      // Constructors
      //

      BandSVDiv(const GenBandMatrix<T>& A, bool StoreU, bool StoreV);
      ~BandSVDiv();

      //
      // Div, DivEq
      //

      template <class T1> void DoLDivEq2(const MatrixView<T1>& m) const;
      template <class T1> inline void DoLDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(Uisset() && Visset());
	DoLDivEq2(m);
      } 

      template <class T1> void DoRDivEq2(const MatrixView<T1>& m) const;
      template <class T1> inline void DoRDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(Uisset() && Visset());
	DoRDivEq2(m);
      } 

      template <class T1, class T2> void DoLDiv2(const GenMatrix<T1>& m,
	    const MatrixView<T2>& x) const;
      template <class T1, class T2> inline void DoLDiv(const GenMatrix<T1>& m,
	    const MatrixView<T2>& x) const 
      {
	TMVAssert(Uisset() && Visset());
	DoLDiv2(m,x);
      } 

      template <class T1, class T2> void DoRDiv2(const GenMatrix<T1>& m,
	    const MatrixView<T2>& x) const;
      template <class T1, class T2> inline void DoRDiv(const GenMatrix<T1>& m,
	    const MatrixView<T2>& x) const 
      {
	TMVAssert(Uisset() && Visset());
	DoRDiv2(m,x);
      } 

#include "TMV_AuxAllDiv.h"

      //
      // Determinant, Inverse
      //

      T Det() const;
      template <class T1> void DoInverse2(const MatrixView<T1>& minv) const;
      template <class T1> inline void DoInverse(
	  const MatrixView<T1>& minv) const
      {
	TMVAssert(Uisset() && Visset());
	DoInverse2(minv);
      }
      void DoInverseATA2(const MatrixView<T>& minv) const;
      inline void DoInverseATA(const MatrixView<T>& minv) const
      {
	TMVAssert(Visset());
	DoInverseATA2(minv);
      }

      bool Singular() const;
      RealType(T) Norm2() const;
      RealType(T) Condition() const;

      //
      // Determine which (if any) S values to zero out
      //

      void Thresh(RealType(T) toler, std::ostream* debugout=0) const;
      void Top(size_t neigen, std::ostream* debugout=0) const;
      size_t GetKMax() const;

      //
      // Access Decomposition
      //

      ConstMatrixView<T> GetU() const;
      ConstVectorView<RealType(T)> GetS() const;
      ConstMatrixView<T> GetV() const;

      std::string Type() const;
      DivType GetDivType() const;
      bool CheckDecomp(const BaseMatrix<T>& m, std::ostream* fout) const;

    protected:

      struct BandSVDiv_Impl;
      BandSVDiv_Impl* pimpl;

      size_t colsize() const;
      size_t rowsize() const;
      bool Visset() const;
      bool Uisset() const;

  }; // BandSVDiv

} // namespace tmv;

#endif
