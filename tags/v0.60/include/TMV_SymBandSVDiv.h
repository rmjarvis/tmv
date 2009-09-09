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
#ifndef TMV_SymBandSVDiv_H
#define TMV_SymBandSVDiv_H

#include "TMV_SymDivider.h"

namespace tmv {

  template <class T> void HermBandSV_Decompose(
      const GenSymBandMatrix<T>& A,
      const MatrixView<T>* U, const VectorView<RealType(T)>& S, T& det);

  template <class T> void SymBandSV_Decompose(
      const GenSymBandMatrix<T>& A,
      const MatrixView<T>* U, const VectorView<RealType(T)>& S, 
      const MatrixView<T>* V, T& det);

  template <class T> class HermBandSVDiv : 
    public SymDivider<T> 
  {

    public :

      //
      // Constructors
      //

      HermBandSVDiv(const GenSymBandMatrix<T>& A, bool StoreU);
      ~HermBandSVDiv();

      //
      // Div, DivEq
      //

      template <class T1> void DoLDivEq2(const MatrixView<T1>& m) const;
      template <class T1> void DoLDivEq(const MatrixView<T1>& m) const
      {
	TMVAssert(Uisset());
	DoLDivEq2(m);
      }
      template <class T1> void DoRDivEq2(const MatrixView<T1>& m) const;
      template <class T1> void DoRDivEq(const MatrixView<T1>& m) const
      {
	TMVAssert(Uisset());
	DoRDivEq2(m);
      }
      template <class T1, class T2> void DoLDiv2(const GenMatrix<T1>& m,
	  const MatrixView<T2>& x) const;
      template <class T1, class T2> void DoLDiv(const GenMatrix<T1>& m,
	  const MatrixView<T2>& x) const
      {
	TMVAssert(Uisset());
	DoLDiv2(m,x);
      }
      template <class T1, class T2> void DoRDiv2(const GenMatrix<T1>& m, 
	  const MatrixView<T2>& x) const;
      template <class T1, class T2> void DoRDiv(const GenMatrix<T1>& m, 
	  const MatrixView<T2>& x) const
      {
	TMVAssert(Uisset());
	DoRDiv2(m,x);
      }

      //
      // Determinant, Inverse
      //

      T Det() const;
      template <class T1> void DoInverse2(const MatrixView<T1>& minv) const;
      template <class T1> void DoInverse(const MatrixView<T1>& minv) const
      {
	TMVAssert(Uisset());
	DoInverse2(minv);
      }
      template <class T1> void DoInverse(const SymMatrixView<T1>& sinv) const;
      inline void Inverse(const SymMatrixView<RealType(T)>& sinv) const
      {
	TMVAssert(Uisset());
	TMVAssert(IsReal(T()));
	TMVAssert(sinv.size() == colsize());
	DoInverse(sinv);
      }
      inline void Inverse(const SymMatrixView<ComplexType(T)>& sinv) const
      {
	TMVAssert(Uisset());
	TMVAssert(sinv.size() == colsize());
	TMVAssert(sinv.isherm());
	DoInverse(sinv);
      }
      void DoInverseATA2(const MatrixView<T>& minv) const;
      inline void DoInverseATA(const MatrixView<T>& minv) const
      { 
	TMVAssert(Uisset());
	DoInverseATA2(minv);
      }
      bool Singular() const;
      RealType(T) Norm2() const;
      RealType(T) Condition() const;

#include "TMV_AuxAllDiv.h"

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

      struct HermBandSVDiv_Impl;
      HermBandSVDiv_Impl* pimpl;

      size_t colsize() const;
      size_t rowsize() const;
      bool Uisset() const;

  }; // HermSVDiv

  template <class T> class SymBandSVDiv : 
    public SymDivider<T> 
  {

    public :

      //
      // Constructors
      //

      SymBandSVDiv(const GenSymBandMatrix<T>& A, bool StoreU, bool StoreV);
      ~SymBandSVDiv();

      //
      // Div, DivEq
      //

      template <class T1> void DoLDivEq2(const MatrixView<T1>& m) const;
      template <class T1> void DoLDivEq(const MatrixView<T1>& m) const
      {
	TMVAssert(Uisset() && Visset());
	DoLDivEq2(m);
      }
      template <class T1> void DoRDivEq2(const MatrixView<T1>& m) const;
      template <class T1> void DoRDivEq(const MatrixView<T1>& m) const
      {
	TMVAssert(Uisset() && Visset());
	DoRDivEq2(m);
      }
      template <class T1, class T2> void DoLDiv2(
	  const GenMatrix<T1>& m, const MatrixView<T2>& x) const;
      template <class T1, class T2> void DoLDiv(
	  const GenMatrix<T1>& m, const MatrixView<T2>& x) const
      {
	TMVAssert(Uisset() && Visset());
	DoLDiv2(m,x);
      }
      template <class T1, class T2> void DoRDiv2(
	  const GenMatrix<T1>& m, const MatrixView<T2>& x) const;
      template <class T1, class T2> void DoRDiv(
	  const GenMatrix<T1>& m, const MatrixView<T2>& x) const
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
      template <class T1> void DoInverse(const SymMatrixView<T1>& sinv) const;
      inline void Inverse(const SymMatrixView<RealType(T)>& sinv) const
      {
	TMVAssert(Uisset() && Visset());
	TMVAssert(IsReal(T()));
	TMVAssert(sinv.size() == colsize());
	DoInverse(sinv);
      }
      inline void Inverse(const SymMatrixView<ComplexType(T)>& sinv) const
      {
	TMVAssert(Uisset() && Visset());
	TMVAssert(sinv.size() == colsize());
	TMVAssert(sinv.issym());
	DoInverse(sinv);
      }
      void DoInverseATA2(const MatrixView<T>& minv) const;
      inline void DoInverseATA(const MatrixView<T>& minv) const
      {
	TMVAssert(Visset());
	DoInverse2(minv);
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

      struct SymBandSVDiv_Impl;
      SymBandSVDiv_Impl* pimpl;

      size_t colsize() const;
      size_t rowsize() const;
      bool Uisset() const;
      bool Visset() const;

  }; // SymBandSVDiv

} // namespace tmv;

#endif
