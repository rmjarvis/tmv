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
// This file defines the TMV BaseMatrix class.
//
// This base class defines some of the things that all 
// matrices need to be able to do, as well as some of the
// arithmetic operations (those that return a Vector).
// This should be used as the base class for generic
// matrices as well as any special ones (eg. sparse,
// symmetric, etc.)
//
//

#ifndef TMV_BaseMatrix_H
#define TMV_BaseMatrix_H

#include "TMV_Base.h"
#include "TMV_BaseVector.h"

namespace tmv {

  template <class T> class BaseMatrix;
  template <class T> class GenMatrix;
  template <class T, IndexStyle I=CStyle> class ConstMatrixView;
  template <class T, IndexStyle I=CStyle> class MatrixView;
  template <class T, StorageType S=ColMajor, IndexStyle I=CStyle> class Matrix;
  template <class T> class Divider;

  template <class T> struct AssignableToMatrix
  {
    virtual size_t colsize() const = 0;
    virtual size_t rowsize() const = 0;
    inline size_t ncols() const 
    { return rowsize(); }
    inline size_t nrows() const 
    { return colsize(); }
    inline bool IsSquare() const 
    { return colsize() == rowsize(); }

    virtual void AssignToM(const MatrixView<RealType(T)>& m) const = 0; 
    virtual void AssignToM(const MatrixView<ComplexType(T)>& m) const = 0; 

    virtual inline ~AssignableToMatrix() {}
  };

  template <class T> class BaseMatrix :
    virtual public AssignableToMatrix<T>
  {
    public :
      //
      // Access Functions
      //

      using AssignableToMatrix<T>::colsize;
      using AssignableToMatrix<T>::rowsize;

      //
      // Functions of Matrix
      //

      virtual T Det() const = 0;
      virtual RealType(T) LogDet(T* sign=0) const = 0;
      virtual T Trace() const = 0;

      virtual RealType(T) Norm() const  = 0;
      virtual RealType(T) NormSq(RealType(T) scale = RealType(T)(1)) const = 0;
      virtual RealType(T) NormF() const  = 0;
      virtual RealType(T) Norm1() const = 0;
      virtual RealType(T) Norm2() const  = 0;
      virtual RealType(T) NormInf() const = 0;
      virtual RealType(T) MaxAbsElement() const = 0;

      virtual void Inverse(const MatrixView<T>& minv) const = 0;
      virtual void InverseATA(const MatrixView<T>& ata) const = 0;
      virtual bool Singular() const = 0;
      virtual RealType(T) Condition() const = 0;

      virtual auto_ptr<BaseMatrix<T> > NewCopy() const = 0;
      virtual auto_ptr<BaseMatrix<T> > NewView() const = 0;
      virtual auto_ptr<BaseMatrix<T> > NewTranspose() const = 0;
      virtual auto_ptr<BaseMatrix<T> > NewConjugate() const = 0;
      virtual auto_ptr<BaseMatrix<T> > NewAdjoint() const = 0;
      virtual auto_ptr<BaseMatrix<T> > NewInverse() const = 0;

      // 
      // I/O: Write
      //

      virtual void Write(std::ostream& os) const = 0;
      virtual void Write(std::ostream& os, RealType(T) thresh) const = 0;

      virtual inline ~BaseMatrix() {}
  }; // BaseMatrix

  template <class T> class DivHelper :
    virtual public AssignableToMatrix<T>
  {

    public:

      //
      // Constructors
      //

      inline DivHelper() : pdiv(0) {}
      inline DivHelper(const BaseMatrix<T>&) : pdiv(0) {}
      // Cannot do this inline, since need to delete pdiv, 
      // and I only define DivImpl in BaseMatrix.cpp.
      virtual ~DivHelper();

      using AssignableToMatrix<T>::colsize;
      using AssignableToMatrix<T>::rowsize;

      T Det() const 
      {
	TMVAssert(rowsize() == colsize());
	return DoDet(); 
      }

      RealType(T) LogDet(T* sign) const 
      {
	TMVAssert(rowsize() == colsize());
	return DoLogDet(sign); 
      }

      void Inverse(const MatrixView<T>& minv) const
      {
	TMVAssert(minv.colsize() == rowsize());
	TMVAssert(minv.rowsize() == colsize());
	DoInverse(minv);
      }

      template <class T1> inline void Inverse(const MatrixView<T1>& minv) const
      {
	TMVAssert(minv.colsize() == rowsize());
	TMVAssert(minv.rowsize() == colsize());
	DoInverse(minv);
      }

      template <class T1, StorageType S, IndexStyle I> inline void Inverse(
	  Matrix<T1,S,I>& minv) const
      {
	TMVAssert(minv.colsize() == rowsize());
	TMVAssert(minv.rowsize() == colsize());
	DoInverse(minv.View());
      }

      inline void InverseATA(const MatrixView<T>& ata) const
      { 
	TMVAssert(ata.colsize() == 
	    rowsize() < colsize() ? rowsize() : colsize());
	TMVAssert(ata.rowsize() == 
	    rowsize() < colsize() ? rowsize() : colsize());
	DoInverseATA(ata);
      }

      template <StorageType S, IndexStyle I> inline void InverseATA(
	  Matrix<T,S,I>& ata) const
      { 
	TMVAssert(ata.colsize() == 
	    rowsize() < colsize() ? rowsize() : colsize());
	TMVAssert(ata.rowsize() == 
	    rowsize() < colsize() ? rowsize() : colsize());
	DoInverseATA(ata.View());
      }

      inline bool Singular() const
      { return DoSingular(); }

      inline RealType(T) Norm2() const
      {
	TMVAssert(DivIsSet() && GetDivType() == SV);
	return DoNorm2(); 
      }

      inline RealType(T) Condition() const
      {
	TMVAssert(DivIsSet() && GetDivType() == SV);
	return DoCondition(); 
      }

      // m^-1 * v -> v
      template <class T1> inline void LDivEq(const VectorView<T1>& v) const 
      {
	TMVAssert(colsize() == rowsize());
	TMVAssert(colsize() == v.size());
	DoLDivEq(v);
      }

      template <class T1> inline void LDivEq(const MatrixView<T1>& m) const 
      { 
	TMVAssert(colsize() == rowsize());
	TMVAssert(colsize() == m.colsize());
	DoLDivEq(m);
      }

      // v * m^-1 -> v
      template <class T1> inline void RDivEq(const VectorView<T1>& v) const 
      { 
	TMVAssert(colsize() == rowsize());
	TMVAssert(colsize() == v.size());
	DoRDivEq(v);
      }

      template <class T1> inline void RDivEq(const MatrixView<T1>& m) const 
      { 
	TMVAssert(colsize() == rowsize());
	TMVAssert(colsize() == m.rowsize());
	DoRDivEq(m);
      }

      // m^-1 * v1 -> v0
      template <class T1, class T0> inline void LDiv(
	  const GenVector<T1>& v1, const VectorView<T0>& v0) const
      { 
	TMVAssert(rowsize() == v0.size());
	TMVAssert(colsize() == v1.size());
	DoLDiv(v1,v0);
      }
     
      template <class T1, class T0> inline void LDiv(
	  const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
      { 
	TMVAssert(rowsize() == m0.colsize());
	TMVAssert(colsize() == m1.colsize());
	TMVAssert(m1.rowsize() == m0.rowsize());
	DoLDiv(m1,m0);
      }

      // v1 * m^-1 -> v0
      template <class T1, class T0> inline void RDiv(
	  const GenVector<T1>& v1, const VectorView<T0>& v0) const
      { 
	TMVAssert(rowsize() == v1.size());
	TMVAssert(colsize() == v0.size());
	DoRDiv(v1,v0);
      }

      template <class T1, class T0> inline void RDiv(
	  const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
      { 
	TMVAssert(rowsize() == m1.rowsize());
	TMVAssert(colsize() == m0.rowsize());
	TMVAssert(m1.colsize() == m0.colsize());
	DoRDiv(m1,m0);
      }

      //
      // Division Control
      //

      void DivideInPlace() const;
      void SaveDiv() const;
      void DivideUsing(DivType dt) const;
      void SetDiv() const;
      bool DivIsSet() const;
      void UnSetDiv() const;
      void ReSetDiv() const;

      bool CheckDecomp(std::ostream* fout=0) const;
      bool CheckDecomp(const BaseMatrix<T>& m2, std::ostream* fout=0) const;

    protected :

      struct DivImpl;

      mutable DivImpl* pdiv;

      const Divider<T>* GetDiv() const;
      void SetDiv(Divider<T>*) const;
      DivType GetDivType() const;
      bool IsDivInPlace() const;
      void DoneDiv() const;

      // This is why the divider stuff is implemented using private
      // inheritance.  NewDivider needs to be defined in the 
      // derived class.
      virtual void NewDivider() const = 0;
      virtual const BaseMatrix<T>& GetMatrix() const = 0;

    private :

      inline DivHelper(const DivHelper<T>&) : pdiv(0) 
      { TMVAssert(FALSE); }
      inline DivHelper<T>& operator=(const DivHelper<T>&)
      { TMVAssert(FALSE); return *this; }

      void SetupDiv() const;

      T DoDet() const;
      RealType(T) DoLogDet(T* sign) const;
      template <class T1> void DoInverse(const MatrixView<T1>& minv) const;
      void DoInverseATA(const MatrixView<T>& minv) const;
      bool DoSingular() const;
      RealType(T) DoNorm2() const;
      RealType(T) DoCondition() const;
      void DoWrite(std::ostream& os) const;
      void DoWrite(std::ostream& os, RealType(T) thresh) const;
      template <class T1> void DoLDivEq(const VectorView<T1>& v) const;
      template <class T1> void DoLDivEq(const MatrixView<T1>& m) const;
      template <class T1> void DoRDivEq(const VectorView<T1>& v) const;
      template <class T1> void DoRDivEq(const MatrixView<T1>& m) const;
      template <class T1, class T0> void DoLDiv(
	  const GenVector<T1>& v1, const VectorView<T0>& v0) const;
      template <class T1, class T0> void DoLDiv(
	  const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const;
      template <class T1, class T0> void DoRDiv(
	  const GenVector<T1>& v1, const VectorView<T0>& v0) const;
      template <class T1, class T0> void DoRDiv(
	  const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const;

  }; // DivHelper

  //
  // Functions of Matrices:
  //

  template <class T> inline T Det(const BaseMatrix<T>& m)
  { return m.Det(); }

  template <class T> inline RealType(T) LogDet(const BaseMatrix<T>& m)
  { return m.LogDet(); }

  template <class T> inline T Trace(const BaseMatrix<T>& m)
  { return m.Trace(); }

  template <class T> inline RealType(T) Norm(const BaseMatrix<T>& m)
  { return m.Norm(); }

  template <class T> inline RealType(T) NormSq(const BaseMatrix<T>& m)
  { return m.NormSq(); }
  
  template <class T> inline RealType(T) NormF(const BaseMatrix<T>& m)
  { return m.NormF(); }

  template <class T> inline RealType(T) Norm1(const BaseMatrix<T>& m)
  { return m.Norm1(); }

  template <class T> inline RealType(T) Norm2(const BaseMatrix<T>& m)
  { return m.Norm2(); }

  template <class T> inline RealType(T) NormInf(const BaseMatrix<T>& m)
  { return m.NormInf(); }

  template <class T> inline RealType(T) MaxAbsElement(const BaseMatrix<T>& m)
  { return m.MaxAbsElement(); }

  //
  // I/O
  //

  template <class T> inline std::ostream& operator<<(std::ostream& os, 
      const BaseMatrix<T>& m)
  { m.Write(os); return os; }

#ifdef XDEBUG
  inline std::string Text(StorageType s)
  {
    return s == RowMajor ? "RowMajor" :
      s == ColMajor ? "ColMajor" :
      s == DiagMajor ? "DiagMajor" :
      s == NoMajor ? "NoMajor" : "Unknown";
  }

  template <class T, StorageType S, IndexStyle I> inline std::string Type(
      const Matrix<T,S,I>& )
  {
    return std::string("Matrix<")+Type(T())+","+Text(S)+","+Text(I)+">";
  }
  template <class T> inline std::string Type(const GenMatrix<T>& m)
  {
    return std::string("GenMatrix<")+Type(T())+","+Text(m.stor())+
      ","+Text(m.ct())+">";
  }
  template <class T, IndexStyle I> inline std::string Type(
      const ConstMatrixView<T,I>& m)
  {
    return std::string("ConstMatrixView<")+Type(T())+","+Text(m.stor())+","
      +Text(I)+","+Text(m.ct())+">";
  }
  template <class T, IndexStyle I> inline std::string Type(
      const MatrixView<T,I>& m)
  {
    return std::string("MatrixView<")+Type(T())+","+Text(m.stor())+","
      +Text(I)+","+Text(m.ct())+">";
  }
#endif

} // namespace tmv

#endif
