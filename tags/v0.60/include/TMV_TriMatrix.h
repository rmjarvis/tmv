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
//
// This file defines the TMV TriMatrix class.
//
// Constructors:
//
//    There are two TriMatrix classes: UpperTriMatrix<T> and LowerTriMatrix<T>
//    For these notes, I will just write TriMatrix, but for all uses,
//    you need to write "Upper" or "Lower" before the "Tri".
//
//    In addition to the type template parameter (T), TriMatrixes have two
//    additional template parameters:
//        DiagType dt = UnitDiag || NonUnitDiag 
//        StorageType stor = RowMajor || ColMajor
//
//        They both have default values, so you can omit both, or
//        just stor.  The default values are: {NonUnitDiag, RowMajor}
//
//        If dt is UnitDiag, then the diagonal elements are not
//        actually stored or referenced.  The are all taken to be = 1.
//
//        The storage follows the same meaning as for regular Matrices.
//
//    TriMatrix<T,dt,stor>(size_t n)
//        Makes a Triangular Matrix with column size = row size = n
//        with _uninitialized_ values.
//
//    TriMatrix<T,dt,stor>(size_t n, T x)
//        Makes a Triangular Matrix with column size = row size = n
//        with all values = x
//
//    TriMatrix<T,dt,stor>(const Matrix<T>& m)
//    TriMatrix<T,dt,stor>(const TriMatrix<T>& m)
//        Makes a TriMatrix which copies the corresponding elements of m.
//
//    ConstUpperTriMatrixView<T> UpperTriMatrixViewOf(const Matrix<T>& m, dt)
//    ConstLowerTriMatrixView<T> LowerTriMatrixViewOf(const Matrix<T>& m, dt)
//    ConstUpperTriMatrixView<T> UpperTriMatrixViewOf(const TriMatrix<T>& m, dt)
//    ConstLowerTriMatrixView<T> LowerTriMatrixViewOf(const TriMatrix<T>& m, dt)
//        Makes a constant TriMatrix view of the corresponding part of m.
//        While this view cannot be modified, changing the original m
//        will cause corresponding changes in this view of m.
//        The last two versions allows you to re-view a NonUnitDiag TriMatrix
//        as a UnitDiag TriMatrix.
//
//    UpperTriMatrixView<T> UpperTriMatrixViewOf(Matrix<T>& m, dt)
//    LowerTriMatrixView<T> LowerTriMatrixViewOf(Matrix<T>& m, dt)
//    UpperTriMatrixView<T> UpperTriMatrixViewOf(UpperTriMatrix<T>& m, dt)
//    LowerTriMatrixView<T> LowerTriMatrixViewOf(LowerTriMatrix<T>& m, dt)
//        Makes a modifiable TriMatrix view of the corresponding part of m.
//        Modifying this matrix will change the corresponding elements in
//        the original Matrix.
//
//    ConstUpperTriMatrixView<T> UpperTriMatrixViewOf(const T* m, size, 
//        dt, stor)
//    UpperTriMatrixView<T> UpperTriMatrixViewOf(T* m, size, dt, stor)
//    ConstLowerTriMatrixView<T> LowerTriMatrixViewOf(const T* m, size, 
//        dt, stor)
//    LowerTriMatrixView<T> LowerTriMatrixViewOf(T* m, size, dt, stor)
//
// Access Functions
//
//    size_t colsize() const
//    size_t rowsize() const
//    size_t size() const
//        Return the dimensions of the TriMatrix
//
//    T& operator()(size_t i, size_t j)
//    T operator()(size_t i, size_t j) const
//        Return the (i,j) element of the TriMatrix
//
//    Vector& row(size_t i, size_t j1, size_t j2)
//        Return a portion of the ith row 
//        This range must be a valid range for the requested row.
//
//    Vector& col(size_t j, size_t i1, size_t i2)
//        Return a portion of the jth column
//        This range must be a valid range for the requested column.
//
//    Vector& diag()
//        Return the main diagonal
//        The TriMatrix must be NonUnitDiag.
//
//    Vector& diag(int i, size_t j1, size_t j2)
//    Vector& diag(int i)
//        Return the super- or sub-diagonal i
//        If i > 0 return the super diagonal starting at m_0i
//        If i < 0 return the sub diagonal starting at m_|i|0
//        If j1,j2 are given, it returns the diagonal SubVector 
//        either from m_j1,i+j1 to m_j2,i+j2 (for i>0) 
//        or from m_|i|+j1,j1 to m_|i|+j2,j2 (for i<0)
//        i>0 will give an error for a LowerTriMatrix
//        i<0 will give an error for an UpperTriMatrix
//        i=0 will give an error for a UnitDiag TriMatrix
//
// Modifying Functions:
//
//    Zero()
//    SetAllTo(T x)
//    ConjugateSelf()
//    SetToIdentity(x = 1)
//    void Swap(TriMatrix& m1, TriMatrix& m2)
//        The TriMatrices must be the same size and shape (Upper or Lower).
//
// Views of a TriMatrix:
//
//    SubMatrix(int i1, int i2, int j1, int j2, int istep=1, int jstep=1)
//        This member function will return a submatrix using rows i1 to i2
//        and columns j1 to j2 which refers
//        to the same physical elements as the original.
//        The submatrix must be completely contained within the TriMatrix.
//
//    SubVector(int i, int j, int istep, int jstep, int size)
//        Returns a SubVector which starts at position (i,j) in the 
//        matrix, moves in the directions (istep,jstep) and has a length
//        of size.
//
//    SubTriMatrix(int i1, int i2, int istep)
//        Returns the TriMatrix which runs from i1 to i2 along the diagonal
//        (not including i2) with an optional step, and includes the 
//        off diagonal in the same rows/cols.
//
//        For example, with an UpperTriMatrix of size 10, the x's below
//        are the original data, the O's are the SubTriMatrix returned
//        with the command SubTriMatrix(3,11,2), and the #'s are the 
//        SubTriMatrix returned with SubTriMatrix(0,3)
//
//        ###xxxxxxx
//         ##xxxxxxx
//          #xxxxxxx
//           OxOxOxO
//            xxxxxx
//             OxOxO
//              xxxx
//               OxO
//                xx
//                 O
//
//    OffDiag()
//        Returns the (NonUnitDiag) TriMatrix of all the off-diagonal
//        elements of a TriMatrix.
//
//    MakeUnitDiag()
//        Recasts the TriMatrix to have unit diagonal
//
//    Transpose(m)
//    Adjoint(m)
//        Note that the Transpose or Adjoint of an UpperTriMatrix returns 
//        a view which is a LowerTriMatrix, and vice versa.
//    Conjugate(m)
//
//    Real(), Imag()
//        For a complex TriMatrix, returns the real or imaginary part
//        as a real TriMatrix.
//
//
// Functions of Matrixs:
//
//    Det(m)
//    Trace(m)
//    Norm(m) or NormF(m)
//    NormSq(m)
//    Norm1(m) 
//    Norm2(m) 
//    NormInf(m) 
//
//    m.Inverse()
//    Inverse(m)
//    m.Inverse(minv) // Takes either a TriMatrix or Matrix argument
//    m.InverseATA(invata)
//
//    m.NewTranspose()
//    m.NewConjugate()
//    m.NewAdjoint()
//    m.NewInverse()
//    m.NewCopy()
//
//
// I/O: 
//
//    os << m 
//        Writes m to ostream os in the usual Matrix format
//
//    m.WriteCompact(os)
//        Writes m to ostream os in the following compact format:
//        For an UpperTriMatrix:
//          size 
//          ( m(0,0) m(0,1) ... m(0,size) )
//          ( m(1,1) .. m(1,size) )
//          ...
//          ( m(size,size) )
//
//        For a LowerTriMatrix:
//          size 
//          ( m(0,0) )
//          ( m(1,0) m(1,1) )
//          ...
//          ( m(size,0) ... m(size,size) )
//
//    is >> m
//        Reads m from istream is in the compact format
//        m must already be the correct size for this to work.
//
//    is >> mptr
//        If you do not know the size of the TriMatrix to be read, you can
//        use this form where mptr is an auto_ptr to an undefined TriMatrix.
//        (Note: if the DiagType for the TriMatrix is UnitDiag, then
//        all of the diagonals read in must be = 1.)
//
//
// Division Control Functions:
//
//    Most of the point of using TriMatrixes is that they are easy
//    to divide using either forward substitution or back substitution.
//    Therefore, the only division available for TriMatrixes is 
//    this variety.  To do something else (like SVD), you need to 
//    copy it to a regular matrix.
//


#ifndef TMV_TriMatrix_H
#define TMV_TriMatrix_H

#include "TMV_BaseTriMatrix.h"
#include "TMV_Matrix.h"
#include "TMV_DiagMatrix.h"

namespace tmv {

  template <class T> class GenUpperTriMatrix : 
    virtual public AssignableToUpperTriMatrix<T>,
    public BaseMatrix<T>
  {

    public:

      //
      // Constructors
      //

      inline GenUpperTriMatrix() {}
      inline GenUpperTriMatrix(const GenUpperTriMatrix<T>&) {}
      virtual inline ~GenUpperTriMatrix() {}

      //
      // Access Functions
      //

      using AssignableToUpperTriMatrix<T>::size;
      using AssignableToUpperTriMatrix<T>::dt;
      inline size_t colsize() const { return size(); }
      inline size_t rowsize() const { return size(); }

      inline T operator()(size_t i, size_t j) const
      {
	TMVAssert(i<size());
	TMVAssert(j<size());
	if (i>j) return T(0);
	else if (isunit() && i==j) return T(1);
	else {
	  TMVAssert(okij(i,j));
	  return cref(i,j);
	}
      }

      inline ConstVectorView<T> row(size_t i, size_t j1, size_t j2) const 
      { 
	TMVAssert(i<size());
	TMVAssert(j1<=j2);
	TMVAssert(j2<=size());
	TMVAssert(j1==j2 || okij(i,j1));
	return ConstVectorView<T>(cptr()+i*stepi()+j1*stepj(),j2-j1,stepj(),
	    ct()); 
      }

      inline ConstVectorView<T> col(size_t j, size_t i1, size_t i2) const
      {
	TMVAssert(j<size());
	TMVAssert(i1<=i2);
	TMVAssert(i2<=size());
	TMVAssert(i1==i2 || okij(i2-1,j));
	return ConstVectorView<T>(cptr()+i1*stepi()+j*stepj(),i2-i1,stepi(),
	    ct()); 
      }

      inline ConstVectorView<T> diag() const
      {
	TMVAssert(!isunit());
	return ConstVectorView<T>(cptr(),size(),stepi()+stepj(),ct()); 
      }

      inline ConstVectorView<T> diag(int i) const
      {
	TMVAssert(i<=int(size())); 
	TMVAssert(isunit() ? i>0 : i>=0);
	return ConstVectorView<T>(cptr()+i*stepj(),size()-i,stepi()+stepj(),
	    ct()); 
      }

      inline ConstVectorView<T> diag(int i, size_t j1, size_t j2) const
      {
	TMVAssert(isunit() ? i>0 : i>=0);
	TMVAssert(i<=int(size())); 
	TMVAssert(j1 <= j2);
	TMVAssert(j2 <= size()-i);
	const int ds = stepi()+stepj();
	return ConstVectorView<T>(cptr()+i*stepj()+j1*ds,j2-j1,ds,ct());
      }

      template <class T2> inline bool SameAs(const BaseMatrix<T2>& ) const
      { return false; }

      inline bool SameAs(const GenUpperTriMatrix<T>& m2) const
      { 
	return (this == &m2 || (cptr()==m2.cptr() && size()==m2.size() && 
	    dt() == m2.dt() && ct() == m2.ct() &&
	    stepi()==m2.stepi() && stepj()==m2.stepj()));
      }

      inline void AssignToM(const MatrixView<RealType(T)>& m2) const
      {
	TMVAssert(m2.colsize() == size());
	TMVAssert(m2.rowsize() == size());
	TMVAssert(IsReal(T()));
	AssignToU(UpperTriMatrixViewOf(m2,dt()));
	if (isunit()) m2.diag().SetAllTo(RealType(T)(1));
	if (size() > 0) LowerTriMatrixViewOf(m2).OffDiag().Zero();
      }

      inline void AssignToM(const MatrixView<ComplexType(T)>& m2) const
      {
	TMVAssert(m2.colsize() == size());
	TMVAssert(m2.rowsize() == size());
	AssignToU(UpperTriMatrixViewOf(m2,dt()));
	if (isunit()) m2.diag().SetAllTo(T(1));
	if (size() > 0) LowerTriMatrixViewOf(m2).OffDiag().Zero();
      }

      inline void AssignToU(const UpperTriMatrixView<RealType(T)>& m2) const
      { 
	TMVAssert(m2.size() == size());
	TMVAssert(isunit() || !m2.isunit());
	TMVAssert(IsReal(T()));
	if (!SameAs(m2)) Copy(*this,m2); 
      }

      inline void AssignToU(const UpperTriMatrixView<ComplexType(T)>& m2) const
      { 
	TMVAssert(m2.size() == size());
	TMVAssert(isunit() || !m2.isunit());
	if (!SameAs(m2)) Copy(*this,m2); 
      }

      //
      // SubMatrix
      //

      bool OKSubMatrix(
	  int i1, int i2, int j1, int j2, int istep, int jstep) const;

      bool OKSubVector(
	  int i, int j, int istep, int jstep, size_t size) const;

      bool OKSubTriMatrix(int i1, int i2, int istep) const;

      inline ConstMatrixView<T> SubMatrix(
	  int i1, int i2, int j1, int j2) const
      {
	TMVAssert(OKSubMatrix(i1,i2,j1,j2,1,1));
	return ConstMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, stepi(), stepj(), stor(), ct());
      }

      inline ConstMatrixView<T> SubMatrix(
	  int i1, int i2, int j1, int j2, int istep, int jstep) const
      {
	const StorageType newstor =
	  iscm() ? (istep == 1 ? ColMajor : NoMajor) :
	    isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
	TMVAssert(OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	return ConstMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
	    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(), 
	    newstor, ct());
      }

      inline ConstVectorView<T> SubVector(
	  int i, int j, int istep, int jstep, size_t size) const
      {
	TMVAssert(OKSubVector(i,j,istep,jstep,size));
	return ConstVectorView<T>(cptr()+i*stepi()+j*stepj(),size,
	    istep*stepi()+jstep*stepj(),ct());
      }

      inline ConstUpperTriMatrixView<T> SubTriMatrix(int i1, int i2) const
      {
	TMVAssert(OKSubTriMatrix(i1,i2,1));
	return ConstUpperTriMatrixView<T>(cptr()+i1*(stepi()+stepj()),
	    i2-i1,stepi(),stepj(),dt(),stor(),ct());
      }

      inline ConstUpperTriMatrixView<T> SubTriMatrix(
	  int i1, int i2, int istep) const
      {
	TMVAssert(OKSubTriMatrix(i1,i2,istep));
	return ConstUpperTriMatrixView<T>(cptr()+i1*(stepi()+stepj()),
	    (i2-i1)/istep,istep*stepi(),istep*stepj(),dt(),
	    istep==1 ? stor() : NoMajor,ct());
      }

      inline ConstUpperTriMatrixView<T> OffDiag(size_t noff=1) const
      {
	TMVAssert(size() >= noff);
	return ConstUpperTriMatrixView<T>(cptr()+noff*stepj(),size()-noff,
	    stepi(),stepj(),NonUnitDiag,stor(),ct());
      }

      inline ConstUpperTriMatrixView<T> MakeUnitDiag() const
      { 
	return ConstUpperTriMatrixView<T>(cptr(),size(),
	    stepi(),stepj(),UnitDiag,stor(),ct());
      }

      inline ConstUpperTriMatrixView<T> View() const
      { 
	return ConstUpperTriMatrixView<T>(cptr(),size(),
	    stepi(),stepj(),dt(),stor(),ct());
      }

      inline ConstLowerTriMatrixView<T> Transpose() const
      { 
	return ConstLowerTriMatrixView<T>(cptr(),size(),
	    stepj(),stepi(),dt(),TransOf(stor()),ct());
      }

      inline ConstUpperTriMatrixView<T> Conjugate() const
      { 
	return ConstUpperTriMatrixView<T>(cptr(),size(),
	    stepi(),stepj(),dt(),stor(),ConjOf(T,ct()));
      }

      inline ConstLowerTriMatrixView<T> Adjoint() const
      { 
	return ConstLowerTriMatrixView<T>(cptr(),size(),
	    stepj(),stepi(),dt(),TransOf(stor()),ConjOf(T,ct()));
      }

      inline ConstUpperTriMatrixView<RealType(T)> Real() const
      {
	return ConstUpperTriMatrixView<RealType(T)>(
	    reinterpret_cast<const RealType(T)*>(cptr()),size(),
	    IsReal(T()) ? stepi() : 2*stepi(),
	    IsReal(T()) ? stepj() : 2*stepj(),
	    dt(), IsReal(T()) ? stor() : NoMajor,NonConj);
      }

      inline ConstUpperTriMatrixView<RealType(T)> Imag() const
      {
	TMVAssert(IsComplex(T()));
	TMVAssert(!isunit());
	// Since Imag of a UnitDiag TriMatrix has 0's on diagonal.
	return ConstUpperTriMatrixView<RealType(T)>(
	    reinterpret_cast<const RealType(T)*>(cptr())+1,
	    size(),2*stepi(),2*stepj(),dt(),NoMajor,NonConj);
      }

      //
      // Functions of Matrix
      //

      T Det() const;

      inline T Trace() const
      { return isunit() ? T(size()) : diag().SumElements(); }

      inline RealType(T) Norm() const 
      { return SQRT(NormSq()); }

      RealType(T) NormF() const;

      // NormF()^2
      RealType(T) NormSq() const;

      // 1-Norm = max_j (sum_i |a_ij|)
      RealType(T) Norm1() const;

      // 2-Norm needs SVD, so need to do from a regular Matrix
      inline RealType(T) Norm2() const
      { return Matrix<T>(*this).Norm2(); }

      inline RealType(T) Condition() const
      { return Matrix<T>(*this).Condition(); }

      // inf-Norm = max_i (sum_j |a_ij|)
      RealType(T) NormInf() const;

      // = max_i,j (|a_ij|)
      RealType(T) MaxAbsElement() const;

      bool Singular() const { return Det() == T(0); }

      template <class T1> void DoInverse(
	  const UpperTriMatrixView<T1>& minv) const;
      template <class T1> void DoInverse(
	  const MatrixView<T1>& minv) const;
      void DoInverseATA(const MatrixView<T>& ata) const;

      inline void Inverse(const MatrixView<T>& minv) const
      {
	TMVAssert(minv.colsize() == size());
	TMVAssert(minv.rowsize() == size());
	DoInverse(minv);
      }

      template <class T1> inline void Inverse(const MatrixView<T1>& minv) const
      {
	TMVAssert(minv.colsize() == size());
	TMVAssert(minv.rowsize() == size());
	DoInverse(minv);
      }

      template <class T1> inline void Inverse(
	  const UpperTriMatrixView<T1>& minv) const
      {
	TMVAssert(minv.size() == size());
	DoInverse(minv);
      }

      QuotXU<T,T> QInverse() const;
      inline QuotXU<T,T> Inverse() const
      { return QInverse(); }

      inline void InverseATA(const MatrixView<T>& minv) const
      {
	TMVAssert(minv.colsize() == size());
	TMVAssert(minv.rowsize() == size());
	DoInverse(minv);
      }

      template <class T1, DiagType D, StorageType S, IndexStyle I> 
	inline void Inverse(UpperTriMatrix<T1,D,S,I>& minv) const
	{ 
	  TMVAssert(D==NonUnitDiag || isunit());
	  Inverse(minv.View()); 
	}

      template <class T1, StorageType S, IndexStyle I> inline void Inverse(
	  Matrix<T1,S,I>& minv) const
      { Inverse(minv.View()); }

      template <StorageType S, IndexStyle I> inline void InverseATA(
	  Matrix<T,S,I>& minv) const
      { InverseATA(minv.View()); }

      auto_ptr<BaseMatrix<T> > NewCopy() const;
      auto_ptr<BaseMatrix<T> > NewView() const;
      auto_ptr<BaseMatrix<T> > NewTranspose() const;
      auto_ptr<BaseMatrix<T> > NewConjugate() const;
      auto_ptr<BaseMatrix<T> > NewAdjoint() const;
      auto_ptr<BaseMatrix<T> > NewInverse() const;

      //
      // I/O
      //

      void WriteCompact(std::ostream& os) const;
      void Write(std::ostream& os) const;
      void WriteCompact(std::ostream& os, RealType(T) thresh) const;
      void Write(std::ostream& os, RealType(T) thresh) const;

      //
      // Arithmetic Helpers
      //

      template <class T1> void DoLDivEq(const VectorView<T1>& v) const;
      template <class T1, class T0> void DoLDiv(
	  const GenVector<T1>& v1, const VectorView<T0>& v0) const;
      template <class T1> void DoLDivEq(const MatrixView<T1>& m) const;
      template <class T1, class T0> void DoLDiv(
	  const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const;
      template <class T1> void DoLDivEq(const UpperTriMatrixView<T1>& m) const;
      template <class T1, class T0> void DoLDiv(
	  const GenUpperTriMatrix<T1>& m1,
	  const UpperTriMatrixView<T0>& m0) const;

      template <class T1> inline void LDivEq(const VectorView<T1>& v) const
      { 
	TMVAssert(v.size() == size());
	DoLDivEq(v);
      }
      template <class T1, class T0> inline void LDiv(
	  const GenVector<T1>& v1, const VectorView<T0>& v0) const
      { 
	TMVAssert(v0.size() == size());
	TMVAssert(v1.size() == size());
	DoLDiv(v1,v0);
      }
      template <class T1> inline void RDivEq(const VectorView<T1>& v) const
      { Transpose().LDivEq(v); }
      template <class T1, class T0> inline void RDiv(
	  const GenVector<T1>& v1, const VectorView<T0>& v0) const
      { Transpose().LDiv(v1,v0); }

      template <class T1> inline void LDivEq(const MatrixView<T1>& m) const
      { 
	TMVAssert(m.colsize() == size());
	DoLDivEq(m);
      }
      template <class T1, class T0> inline void LDiv(
	  const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
      { 
	TMVAssert(m0.colsize() == size());
	TMVAssert(m1.colsize() == size());
	TMVAssert(m0.rowsize() == m1.rowsize());
	DoLDiv(m1,m0);
      }
      template <class T1> inline void RDivEq(const MatrixView<T1>& m) const
      { Transpose().LDivEq(m.Transpose()); }
      template <class T1, class T0> inline void RDiv(
	  const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
      { Transpose().LDiv(m1.Transpose(),m0.Transpose()); }

      template <class T1> inline void LDivEq(
	  const UpperTriMatrixView<T1>& m) const
      { 
	TMVAssert(m.colsize() == size());
	DoLDivEq(m);
      }
      template <class T1, class T0> inline void LDiv(
	  const GenUpperTriMatrix<T1>& m1,
	  const UpperTriMatrixView<T0>& m0) const
      { 
	TMVAssert(m0.size() == size());
	TMVAssert(m1.size() == size());
	DoLDiv(m1,m0);
      }
      template <class T1> inline void RDivEq(
	  const UpperTriMatrixView<T1>& m) const
      { Transpose().LDivEq(m.Transpose()); }
      template <class T1, class T0> inline void RDiv(
	  const GenUpperTriMatrix<T1>& m1,
	  const UpperTriMatrixView<T0>& m0) const
      { Transpose().LDiv(m1.Transpose(),m0.Transpose()); }

      // For easier compatibility with regular matrices:
      inline void DivideInPlace() const {}
      inline void SaveDiv() const {}
      inline void SetDiv() const {}
      inline void UnSetDiv() const {}
      inline void ReSetDiv() const {}
      inline void DivideUsing(DivType DEBUGPARAM(dt)) const
      { TMVAssert(dt == LU); }

      virtual const T* cptr() const = 0;
      virtual int stepi() const = 0;
      virtual int stepj() const = 0;
      virtual StorageType stor() const = 0;
      virtual ConjType ct() const = 0;
      inline bool isrm() const { return stor() == RowMajor; }
      inline bool iscm() const { return stor() == ColMajor; }
      inline bool isunit() const { return dt() == UnitDiag; }
      inline bool isconj() const
      {
	TMVAssert(IsComplex(T()) || ct()==NonConj);
	return IsComplex(T()) && ct()==Conj;
      }

    protected :

      inline bool okij(size_t i, size_t j) const
      {
	TMVAssert(i < size());
	TMVAssert(j < size());
	if (isunit()) return i<j; else return i<=j;
      }
      virtual T cref(size_t i, size_t j) const;

    private :

      inline void operator=(const GenUpperTriMatrix<T>&) { TMVAssert(FALSE); }

  }; // GenUpperTriMatrix

  template <class T> class GenLowerTriMatrix : 
    virtual public AssignableToLowerTriMatrix<T>,
    public BaseMatrix<T>
  {

    public:

      //
      // Constructors
      //

      inline GenLowerTriMatrix() {}
      inline GenLowerTriMatrix(const GenLowerTriMatrix<T>& ) {}
      virtual inline ~GenLowerTriMatrix() {}

      //
      // Access Functions
      //

      using AssignableToLowerTriMatrix<T>::size;
      using AssignableToLowerTriMatrix<T>::dt;
      inline size_t colsize() const { return size(); }
      inline size_t rowsize() const { return size(); }

      inline T operator()(size_t i, size_t j) const
      {
	TMVAssert(i<size());
	TMVAssert(j<size());
	if (i<j) return T(0);
	else if (isunit() && i==j) return T(1);
	else {
	  TMVAssert(okij(i,j));
	  return cref(i,j);
	}
      }

      inline ConstVectorView<T> row(size_t i, size_t j1, size_t j2) const 
      { 
	TMVAssert(i<size());
	TMVAssert(j1<=j2);
	TMVAssert(j2<=size());
	TMVAssert(j1==j2 || okij(i,j2-1));
	return ConstVectorView<T>(cptr()+i*stepi()+j1*stepj(),
	    j2-j1,stepj(),ct()); 
      }

      inline ConstVectorView<T> col(size_t j, size_t i1, size_t i2) const
      {
	TMVAssert(j<size());
	TMVAssert(i1<=i2);
	TMVAssert(i2<=size());
	TMVAssert(i1==i2 || okij(i1,j));
	return ConstVectorView<T>(cptr()+i1*stepi()+j*stepj(),
	    i2-i1,stepi(),ct()); 
      }

      inline ConstVectorView<T> diag() const
      {
	TMVAssert(!isunit());
	return ConstVectorView<T>(cptr(),size(),stepi()+stepj(),ct()); 
      }

      inline ConstVectorView<T> diag(int i) const
      {
	TMVAssert(-i<=int(size())); 
	TMVAssert(isunit() ? i<0 : i<=0);
	return ConstVectorView<T>(cptr()-i*stepi(),
	    size()-i,stepi()+stepj(),ct());
      }

      inline ConstVectorView<T> diag(int i, size_t j1, size_t j2) const
      {
	TMVAssert(isunit() ? i<0 : i<=0);
	TMVAssert(-i<=int(size())); 
	TMVAssert(j1 <= j2);
	TMVAssert(j2 <= size()+i);
	const int ds = stepi()+stepj();
	return ConstVectorView<T>(cptr()-i*stepi()+j1*ds,j2-j1,ds,ct());
      }

      template <class T2> inline bool SameAs(const BaseMatrix<T2>& ) const
      { return false; }

      inline bool SameAs(const GenLowerTriMatrix<T>& m2) const
      { 
	if (this == &m2) return true;
	else return (cptr()==m2.cptr() && size()==m2.size() && 
	    dt() == m2.dt() && ct() == m2.ct() &&
	    stepi()==m2.stepi() && stepj()==m2.stepj());
      }

      inline void AssignToM(const MatrixView<RealType(T)>& m2) const
      { Transpose().AssignToM(m2.Transpose()); }
      inline void AssignToM(const MatrixView<ComplexType(T)>& m2) const
      { Transpose().AssignToM(m2.Transpose()); }
      inline void AssignToL(const LowerTriMatrixView<RealType(T)>& m2) const
      { Transpose().AssignToU(m2.Transpose()); }
      inline void AssignToL(const LowerTriMatrixView<ComplexType(T)>& m2) const
      { Transpose().AssignToU(m2.Transpose()); }

      //
      // SubMatrix
      //

      inline bool OKSubMatrix(
	  int i1, int i2, int j1, int j2, int istep, int jstep) const
      { return Transpose().OKSubMatrix(j1,j2,i1,i2,jstep,istep); }

      inline bool OKSubVector(
	  int i, int j, int istep, int jstep, size_t size) const
      { return Transpose().OKSubVector(j,i,jstep,istep,size); }

      inline bool OKSubTriMatrix(int i1, int i2, int istep) const
      { return Transpose().OKSubTriMatrix(i1,i2,istep); }

      inline ConstMatrixView<T> SubMatrix(
	  int i1, int i2, int j1, int j2) const
      {
	TMVAssert(OKSubMatrix(i1,i2,j1,j2,1,1));
	return ConstMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, stepi(), stepj(), stor(), ct());
      }

      inline ConstMatrixView<T> SubMatrix(
	  int i1, int i2, int j1, int j2, int istep, int jstep) const
      {
	const StorageType newstor =
	  iscm() ? (istep == 1 ? ColMajor : NoMajor) :
	    isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
	TMVAssert(OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	return ConstMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
	    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(), 
	    newstor,ct());
      }

      inline ConstVectorView<T> SubVector(
	  int i, int j, int istep, int jstep, size_t size) const
      {
	TMVAssert(OKSubVector(i,j,istep,jstep,size));
	return ConstVectorView<T>(cptr()+i*stepi()+j*stepj(),size,
	    istep*stepi()+jstep*stepj(),ct());
      }

      inline ConstLowerTriMatrixView<T> SubTriMatrix(int i1, int i2) const
      {
	TMVAssert(OKSubTriMatrix(i1,i2,1));
	return ConstLowerTriMatrixView<T>(cptr()+i1*(stepi()+stepj()),
	    i2-i1,stepi(),stepj(),dt(),stor(),ct());
      }

      inline ConstLowerTriMatrixView<T> SubTriMatrix(
	  int i1, int i2, int istep) const
      {
	TMVAssert(OKSubTriMatrix(i1,i2,istep));
	return ConstLowerTriMatrixView<T>(cptr()+i1*(stepi()+stepj()),
	    (i2-i1)/istep,istep*stepi(),istep*stepj(),dt(),
	    istep==1 ? stor() : NoMajor,ct());
      }

      inline ConstLowerTriMatrixView<T> OffDiag(size_t noff=1) const
      {
	TMVAssert(size() >= noff);
	return ConstLowerTriMatrixView<T>(cptr()+noff*stepi(),size()-noff,
	    stepi(),stepj(),NonUnitDiag,stor(),ct());
      }

      inline ConstLowerTriMatrixView<T> MakeUnitDiag() const
      { 
	return ConstLowerTriMatrixView<T>(cptr(),size(),
	    stepi(),stepj(),UnitDiag,stor(),ct());
      }

      inline ConstLowerTriMatrixView<T> View() const
      { 
	return ConstLowerTriMatrixView<T>(cptr(),size(),
	    stepi(),stepj(),dt(),stor(),ct());
      }

      inline ConstUpperTriMatrixView<T> Transpose() const
      { 
	return ConstUpperTriMatrixView<T>(cptr(),size(),
	    stepj(),stepi(),dt(),TransOf(stor()),ct());
      }
  
      inline ConstLowerTriMatrixView<T> Conjugate() const
      { 
	return ConstLowerTriMatrixView<T>(cptr(),size(),
	    stepi(),stepj(),dt(),stor(),ConjOf(T,ct()));
      }

      inline ConstUpperTriMatrixView<T> Adjoint() const
      { 
	return ConstUpperTriMatrixView<T>(cptr(),size(),
	    stepj(),stepi(),dt(),TransOf(stor()),ConjOf(T,ct()));
      }

      inline ConstLowerTriMatrixView<RealType(T)> Real() const
      {
	return ConstLowerTriMatrixView<RealType(T)>(
	    reinterpret_cast<const RealType(T)*>(cptr()),size(),
	    IsReal(T()) ? stepi() : 2*stepi(),
	    IsReal(T()) ? stepj() : 2*stepj(),
	    dt(), IsReal(T()) ? stor() : NoMajor,NonConj);
      }

      inline ConstLowerTriMatrixView<RealType(T)> Imag() const
      {
	TMVAssert(IsComplex(T()));
	TMVAssert(!isunit());
	return ConstLowerTriMatrixView<RealType(T)>(
	    reinterpret_cast<const RealType(T)*>(cptr())+1,
	    size(),2*stepi(),2*stepj(),dt(),NoMajor,NonConj);
      }

      //
      // Functions of Matrix
      //

      T Det() const;

      inline T Trace() const
      { return isunit() ? T(size()) : diag().SumElements(); }

      inline RealType(T) Norm() const 
      { return SQRT(NormSq()); }

      inline RealType(T) NormF() const 
      { return SQRT(NormSq()); }

      // NormF()^2
      inline RealType(T) NormSq() const 
      { return Transpose().NormSq(); }

      // 1-Norm = max_j (sum_i |a_ij|)
      inline RealType(T) Norm1() const
      { return Transpose().NormInf(); }

      inline RealType(T) Norm2() const
      { return Matrix<T>(*this).Norm2(); }

      inline RealType(T) Condition() const
      { return Matrix<T>(*this).Condition(); }

      // inf-Norm = max_i (sum_j |a_ij|)
      inline RealType(T) NormInf() const
      { return Transpose().Norm1(); }

      // = max_i,j (|a_ij|)
      inline RealType(T) MaxAbsElement() const
      { return Transpose().MaxAbsElement(); }

      bool Singular() const { return Det() == T(0); }

      QuotXL<T,T> QInverse() const;
      inline QuotXL<T,T> Inverse() const
      { return QInverse(); }

      template <class T1> inline void Inverse(
	  const LowerTriMatrixView<T1>& minv) const
      {
	TMVAssert(minv.size() == size());
	Transpose().Inverse(minv.Transpose());
      }

      inline void Inverse(const MatrixView<T>& minv) const
      {
	TMVAssert(minv.colsize() == size());
	TMVAssert(minv.rowsize() == size());
	Transpose().Inverse(minv.Transpose());
      }

      template <class T1> inline void Inverse(const MatrixView<T1>& minv) const
      {
	TMVAssert(minv.colsize() == size());
	TMVAssert(minv.rowsize() == size());
	Transpose().Inverse(minv.Transpose());
      }

      void DoInverseATA(const MatrixView<T>& minv) const;

      inline void InverseATA(const MatrixView<T>& ata) const
      {
	TMVAssert(ata.colsize() == size());
	TMVAssert(ata.rowsize() == size());
	DoInverseATA(ata);
      }

      template <class T1, DiagType D, StorageType S, IndexStyle I> 
	inline void Inverse(LowerTriMatrix<T1,D,S,I>& minv) const
	{ 
	  TMVAssert(D==NonUnitDiag || isunit());
	  Inverse(minv.View()); 
	}

      template <class T1, StorageType S, IndexStyle I> inline void Inverse(
	  Matrix<T1,S,I>& minv) const
      { Inverse(minv.View()); }

      template <StorageType S, IndexStyle I> inline void InverseATA(
	  Matrix<T,S,I>& minv) const
      { InverseATA(minv.View()); }

      auto_ptr<BaseMatrix<T> > NewCopy() const;
      auto_ptr<BaseMatrix<T> > NewView() const;
      auto_ptr<BaseMatrix<T> > NewTranspose() const ;
      auto_ptr<BaseMatrix<T> > NewConjugate() const;
      auto_ptr<BaseMatrix<T> > NewAdjoint() const;
      auto_ptr<BaseMatrix<T> > NewInverse() const;

      //
      // I/O
      //

      void WriteCompact(std::ostream& os) const;
      void Write(std::ostream& os) const;
      void WriteCompact(std::ostream& os, RealType(T) thresh) const;
      void Write(std::ostream& os, RealType(T) thresh) const;

      //
      // Arithmetic Helpers
      //

      template <class T1> void DoLDivEq(const VectorView<T1>& v) const;
      template <class T1, class T0> void DoLDiv(
	  const GenVector<T1>& v1, const VectorView<T0>& v0) const;
      template <class T1> void DoLDivEq(const MatrixView<T1>& m) const;
      template <class T1, class T0> void DoLDiv(
	  const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const;
      template <class T1> void DoLDivEq(const LowerTriMatrixView<T1>& m) const;
      template <class T1, class T0> void DoLDiv(
	  const GenLowerTriMatrix<T1>& m1,
	  const LowerTriMatrixView<T0>& m0) const;

      template <class T1> inline void LDivEq(const VectorView<T1>& v) const
      { 
	TMVAssert(v.size() == size());
	DoLDivEq(v);
      }
      template <class T1, class T0> inline void LDiv(
	  const GenVector<T1>& v1, const VectorView<T0>& v0) const
      { 
	TMVAssert(v0.size() == size());
	TMVAssert(v1.size() == size());
	DoLDiv(v1,v0);
      }
      template <class T1> inline void RDivEq(const VectorView<T1>& v) const
      { Transpose().LDivEq(v); }
      template <class T1, class T0> inline void RDiv(
	  const GenVector<T1>& v1, const VectorView<T0>& v0) const
      { Transpose().LDiv(v1,v0); }

      template <class T1> inline void LDivEq(const MatrixView<T1>& m) const
      { 
	TMVAssert(m.colsize() == size());
	DoLDivEq(m);
      }
      template <class T1, class T0> inline void LDiv(
	  const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
      { 
	TMVAssert(m0.colsize() == size());
	TMVAssert(m1.colsize() == size());
	TMVAssert(m0.rowsize() == m1.rowsize());
	DoLDiv(m1,m0);
      }
      template <class T1> inline void RDivEq(const MatrixView<T1>& m) const
      { Transpose().LDivEq(m.Transpose()); }
      template <class T1, class T0> inline void RDiv(
	  const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
      { Transpose().LDiv(m1.Transpose(),m0.Transpose()); }

      template <class T1> inline void LDivEq(
	  const LowerTriMatrixView<T1>& m) const
      { 
	TMVAssert(m.colsize() == size());
	DoLDivEq(m);
      }
      template <class T1, class T0> inline void LDiv(
	  const GenLowerTriMatrix<T1>& m1,
	  const LowerTriMatrixView<T0>& m0) const
      { 
	TMVAssert(m0.size() == size());
	TMVAssert(m1.size() == size());
	DoLDiv(m1,m0);
      }
      template <class T1> inline void RDivEq(
	  const LowerTriMatrixView<T1>& m) const
      { Transpose().LDivEq(m.Transpose()); }
      template <class T1, class T0> inline void RDiv(
	  const GenLowerTriMatrix<T1>& m1,
	  const LowerTriMatrixView<T0>& m0) const
      { Transpose().LDiv(m1.Transpose(),m0.Transpose()); }

      // For easier compatibility with regular matrices:
      inline void DivideInPlace() const {}
      inline void SaveDiv() const {}
      inline void SetDiv() const {}
      inline void UnSetDiv() const {}
      inline void ReSetDiv() const {}
      inline void DivideUsing(DivType DEBUGPARAM(dt)) const
      { TMVAssert(dt == LU); }

      virtual const T* cptr() const = 0;
      virtual int stepi() const = 0;
      virtual int stepj() const = 0;
      virtual StorageType stor() const = 0;
      virtual ConjType ct() const = 0;
      inline bool isrm() const { return stor() == RowMajor; }
      inline bool iscm() const { return stor() == ColMajor; }
      inline bool isunit() const { return dt() == UnitDiag; }
      inline bool isconj() const
      {
	TMVAssert(IsComplex(T()) || ct()==NonConj);
	return IsComplex(T()) && ct()==Conj;
      }

    protected :

      inline bool okij(size_t i, size_t j) const
      {
	TMVAssert(i < size());
	TMVAssert(j < size());
	if (isunit()) return i>j; else return i>=j;
      }

      virtual T cref(size_t i, size_t j) const;

    private :

      inline void operator=(const GenLowerTriMatrix<T>&) { TMVAssert(FALSE); }

  }; // GenLowerTriMatrix

  template <class T, IndexStyle I> class ConstUpperTriMatrixView : 
    public GenUpperTriMatrix<T>
  {
    public :

      inline ConstUpperTriMatrixView(const ConstUpperTriMatrixView<T,I>& rhs) :
	itsm(rhs.itsm), itss(rhs.itss), itssi(rhs.itssi), itssj(rhs.itssj),
	itsdiag(rhs.itsdiag), itsstor(rhs.itsstor), itsct(rhs.itsct) {}

      inline ConstUpperTriMatrixView(const GenUpperTriMatrix<T>& rhs) :
	itsm(rhs.cptr()), itss(rhs.size()), 
	itssi(rhs.stepi()), itssj(rhs.stepj()),
	itsdiag(rhs.dt()), itsstor(rhs.stor()), itsct(rhs.ct()) {}

      inline ConstUpperTriMatrixView(
	  const T* _m, size_t _s, int _si, int _sj,
	  DiagType _dt, StorageType _stor, ConjType _ct) : 
	itsm(_m), itss(_s), itssi(_si), itssj(_sj),
	itsdiag(_dt), itsstor(_stor), itsct(_ct)
      { 
	TMVAssert(_stor==RowMajor ? _sj == 1 : _stor==ColMajor ?
	    _si==1 : true);
      }

      virtual inline ~ConstUpperTriMatrixView() {}

      inline size_t size() const { return itss; }
      inline const T* cptr() const { return itsm; }
      inline int stepi() const { return itssi; }
      inline int stepj() const { return itssj; }
      inline StorageType stor() const { return itsstor; }
      inline DiagType dt() const { return itsdiag; }
      inline ConjType ct() const { return itsct; }

    protected :

      const T*const itsm;
      const size_t itss;
      const int itssi;
      const int itssj;
      DiagType itsdiag;
      StorageType itsstor;
      ConjType itsct;

    private :

      inline void operator=(const ConstUpperTriMatrixView<T,I>&) 
      { TMVAssert(FALSE); }

  }; // ConstUpperTriMatrixView

  template <class T, IndexStyle I> class ConstLowerTriMatrixView : 
    public GenLowerTriMatrix<T>
  {
    public :

      inline ConstLowerTriMatrixView(const ConstLowerTriMatrixView<T,I>& rhs) :
	itsm(rhs.itsm), itss(rhs.itss), itssi(rhs.itssi), itssj(rhs.itssj),
	itsdiag(rhs.itsdiag), itsstor(rhs.itsstor), itsct(rhs.itsct) {}

      inline ConstLowerTriMatrixView(const GenLowerTriMatrix<T>& rhs) :
	itsm(rhs.cptr()), itss(rhs.size()), 
	itssi(rhs.stepi()), itssj(rhs.stepj()),
	itsdiag(rhs.dt()), itsstor(rhs.stor()), itsct(rhs.ct()) {}

      inline ConstLowerTriMatrixView(
	  const T* _m, size_t _s, int _si, int _sj,
	  DiagType _dt, StorageType _stor, ConjType _ct) : 
	itsm(_m), itss(_s), itssi(_si), itssj(_sj),
	itsdiag(_dt), itsstor(_stor), itsct(_ct)
      { 
	TMVAssert(_stor==RowMajor ? _sj == 1 : _stor==ColMajor ?
	    _si==1 : true);
      }

      virtual inline ~ConstLowerTriMatrixView() {}

      inline size_t size() const { return itss; }
      inline const T* cptr() const { return itsm; }
      inline int stepi() const { return itssi; }
      inline int stepj() const { return itssj; }
      inline StorageType stor() const { return itsstor; }
      inline DiagType dt() const { return itsdiag; }
      inline ConjType ct() const { return itsct; }

    protected :

      const T*const itsm;
      const size_t itss;
      const int itssi;
      const int itssj;
      DiagType itsdiag;
      StorageType itsstor;
      ConjType itsct;

    private :

      inline void operator=(const ConstLowerTriMatrixView<T,I>&) 
      { TMVAssert(FALSE); }

  }; // ConstLowerTriMatrixView

  template <class T> class ConstUpperTriMatrixView<T,FortranStyle> : 
    public ConstUpperTriMatrixView<T,CStyle>
    {
      public :

	inline ConstUpperTriMatrixView(
	    const ConstUpperTriMatrixView<T,FortranStyle>& rhs) :
	  ConstUpperTriMatrixView<T,CStyle>(rhs) {}

	inline ConstUpperTriMatrixView(
	    const ConstUpperTriMatrixView<T,CStyle>& rhs) :
	  ConstUpperTriMatrixView<T,CStyle>(rhs) {}

	inline ConstUpperTriMatrixView(
	    const GenUpperTriMatrix<T>& rhs) :
	  ConstUpperTriMatrixView<T,CStyle>(rhs) {}

	inline ConstUpperTriMatrixView(
	    const T* _m, size_t _s, int _si, int _sj,
	    DiagType indt, StorageType instor, ConjType inct) : 
	  ConstUpperTriMatrixView<T,CStyle>(_m,_s,_si,_sj,indt,instor,inct) {}

	virtual inline ~ConstUpperTriMatrixView() {}

	//
	// Access Functions
	//

	inline T operator()(size_t i, size_t j) const
	{
	  TMVAssert(i>0 && i<=size());
	  TMVAssert(j>0 && j<=size());
	  if (i>j) return T(0);
	  else if (isunit() && i==j) return T(1);
	  else {
	    TMVAssert(okij(i-1,j-1));
	    return cref(i-1,j-1);
	  }
	}

	inline ConstVectorView<T,FortranStyle> row(
	    size_t i, size_t j1, size_t j2) const 
	{ 
	  TMVAssert(i>0 && i<=size());
	  TMVAssert(j1>0 && j1<=j2 && j2<=size());
	  return GenUpperTriMatrix<T>::row(i-1,j1-1,j2);
	}

	inline ConstVectorView<T,FortranStyle> col(
	    size_t j, size_t i1, size_t i2) const
	{
	  TMVAssert(j>0 && j<=size());
	  TMVAssert(i1>0 && i1<=i2 && i2<=size());
	  return GenUpperTriMatrix<T>::col(j-1,i1-1,i2);
	}

	inline ConstVectorView<T,FortranStyle> diag() const
	{ return GenUpperTriMatrix<T>::diag(); }

	inline ConstVectorView<T,FortranStyle> diag(int i) const
	{ return GenUpperTriMatrix<T>::diag(i); }

	inline ConstVectorView<T,FortranStyle> diag(
	    int i, size_t j1, size_t j2) const
	{
	  TMVAssert(j1>0);
	  return GenUpperTriMatrix<T>::diag(i,j1-1,j2); 
	}

	//
	// SubMatrix
	//

	bool OKSubMatrix(
	    int i1, int i2, int j1, int j2, int istep, int jstep) const;

	bool OKSubVector(
	    int i, int j, int istep, int jstep, size_t size) const;

	bool OKSubTriMatrix(int i1, int i2, int istep) const;

	inline ConstMatrixView<T,FortranStyle> SubMatrix(
	    int i1, int i2, int j1, int j2) const
	{
	  TMVAssert(OKSubMatrix(i1,i2,j1,j2,1,1));
	  return GenUpperTriMatrix<T>::SubMatrix(i1-1,i2,j1-1,j2);
	}

	inline ConstMatrixView<T,FortranStyle> SubMatrix(
	    int i1, int i2, int j1, int j2, int istep, int jstep) const
	{
	  TMVAssert(OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	  return GenUpperTriMatrix<T>::SubMatrix(i1-1,i2-1+istep,
	      j1-1,j2-1+jstep,istep,jstep);
	}

	inline ConstVectorView<T,FortranStyle> SubVector(
	    int i, int j, int istep, int jstep, size_t size) const
	{
	  TMVAssert(OKSubVector(i,j,istep,jstep,size));
	  return GenUpperTriMatrix<T>::SubVector(i-1,j-1,istep,jstep,size);
	}

	inline ConstUpperTriMatrixView<T,FortranStyle> SubTriMatrix(
	    int i1, int i2) const
	{
	  TMVAssert(OKSubTriMatrix(i1,i2,1));
	  return GenUpperTriMatrix<T>::SubTriMatrix(i1-1,i2);
	}

	inline ConstUpperTriMatrixView<T,FortranStyle> SubTriMatrix(
	    int i1, int i2, int istep) const
	{
	  TMVAssert(OKSubTriMatrix(i1,i2,istep));
	  return GenUpperTriMatrix<T>::SubTriMatrix(i1-1,i2-1+istep,istep);
	}

	inline ConstUpperTriMatrixView<T,FortranStyle> OffDiag(
	    size_t noff=1) const
	{ return GenUpperTriMatrix<T>::OffDiag(noff); }

	inline ConstUpperTriMatrixView<T,FortranStyle> MakeUnitDiag() const
	{ return GenUpperTriMatrix<T>::MakeUnitDiag(); }

	inline ConstUpperTriMatrixView<T,FortranStyle> View() const
	{ return GenUpperTriMatrix<T>::View(); }

	inline ConstLowerTriMatrixView<T,FortranStyle> Transpose() const
	{ return GenUpperTriMatrix<T>::Transpose(); }

	inline ConstUpperTriMatrixView<T,FortranStyle> Conjugate() const
	{ return GenUpperTriMatrix<T>::Conjugate(); }

	inline ConstLowerTriMatrixView<T,FortranStyle> Adjoint() const
	{ return GenUpperTriMatrix<T>::Adjoint(); }

	inline ConstUpperTriMatrixView<RealType(T),FortranStyle> Real() const
	{ return GenUpperTriMatrix<T>::Real(); }

	inline ConstUpperTriMatrixView<RealType(T),FortranStyle> Imag() const
	{ return GenUpperTriMatrix<T>::Imag(); }

	using ConstUpperTriMatrixView<T,CStyle>::size;
	using GenUpperTriMatrix<T>::isunit;

      protected :

	using GenUpperTriMatrix<T>::okij;
	using ConstUpperTriMatrixView<T,CStyle>::cref;

      private :

	inline void operator=(const ConstUpperTriMatrixView<T,FortranStyle>&) 
	{ TMVAssert(FALSE); }

    }; // FortranStyle ConstUpperTriMatrixView

  template <class T> class ConstLowerTriMatrixView<T,FortranStyle> : 
    public ConstLowerTriMatrixView<T,CStyle>
    {
      public :

	inline ConstLowerTriMatrixView(
	    const ConstLowerTriMatrixView<T,FortranStyle>& rhs) :
	  ConstLowerTriMatrixView<T,CStyle>(rhs) {}

	inline ConstLowerTriMatrixView(
	    const ConstLowerTriMatrixView<T,CStyle>& rhs) :
	  ConstLowerTriMatrixView<T,CStyle>(rhs) {}

	inline ConstLowerTriMatrixView(
	    const GenLowerTriMatrix<T>& rhs) :
	  ConstLowerTriMatrixView<T,CStyle>(rhs) {}

	inline ConstLowerTriMatrixView(
	    const T* _m, size_t _s, int _si, int _sj,
	    DiagType indt, StorageType instor, ConjType inct) : 
	  ConstLowerTriMatrixView<T,CStyle>(_m,_s,_si,_sj,indt,instor,inct) {}

	virtual inline ~ConstLowerTriMatrixView() {}

	//
	// Access Functions
	//

	inline T operator()(size_t i, size_t j) const
	{
	  TMVAssert(i>0 && i<=size());
	  TMVAssert(j>0 && j<=size());
	  if (i<j) return T(0);
	  else if (isunit() && i==j) return T(1);
	  else {
	    TMVAssert(okij(i-1,j-1));
	    return cref(i-1,j-1);
	  }
	}

	inline ConstVectorView<T,FortranStyle> row(
	    size_t i, size_t j1, size_t j2) const 
	{ 
	  TMVAssert(i>0 && i<=size());
	  TMVAssert(j1>0 && j1<=j2);
	  TMVAssert(j2<=size());
	  return GenLowerTriMatrix<T>::row(i-1,j1-1,j2);
	}

	inline ConstVectorView<T,FortranStyle> col(
	    size_t j, size_t i1, size_t i2) const
	{
	  TMVAssert(j>0 && j<=size());
	  TMVAssert(i1>0 && i1<=i2);
	  TMVAssert(i2<=size());
	  return GenLowerTriMatrix<T>::col(j-1,i1-1,i2);
	}

	inline ConstVectorView<T,FortranStyle> diag() const
	{ return GenLowerTriMatrix<T>::diag(); }

	inline ConstVectorView<T,FortranStyle> diag(int i) const
	{ return GenLowerTriMatrix<T>::diag(i); }

	inline ConstVectorView<T,FortranStyle> diag(
	    int i, size_t j1, size_t j2) const
	{
	  TMVAssert(j1>0);
	  return GenLowerTriMatrix<T>::diag(i,j1-1,j2); 
	}

	//
	// SubMatrix
	//

	inline bool OKSubMatrix(
	    int i1, int i2, int j1, int j2, int istep, int jstep) const
	{ return Transpose().OKSubMatrix(j1,j2,i1,i2,jstep,istep); }

	inline bool OKSubVector(
	    int i, int j, int istep, int jstep, size_t size) const
	{ return Transpose().OKSubVector(j,i,jstep,istep,size); }

	inline bool OKSubTriMatrix(int i1, int i2, int istep) const
	{ return Transpose().OKSubTriMatrix(i1,i2,istep); }

	inline ConstMatrixView<T,FortranStyle> SubMatrix(
	    int i1, int i2, int j1, int j2) const
	{
	  TMVAssert(OKSubMatrix(i1,i2,j1,j2,1,1));
	  return GenLowerTriMatrix<T>::SubMatrix(i1-1,i2,j1-1,j2);
	}

	inline ConstMatrixView<T,FortranStyle> SubMatrix(
	    int i1, int i2, int j1, int j2, int istep, int jstep) const
	{
	  TMVAssert(OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	  return GenLowerTriMatrix<T>::SubMatrix(i1-1,i2-1+istep,
	      j1-1,j2-1+jstep,istep,jstep);
	}

	inline ConstVectorView<T,FortranStyle> SubVector(
	    int i, int j, int istep, int jstep, size_t size) const
	{
	  TMVAssert(OKSubVector(i,j,istep,jstep,size));
	  return GenLowerTriMatrix<T>::SubVector(i-1,j-1,istep,jstep,size);
	}

	inline ConstLowerTriMatrixView<T,FortranStyle> SubTriMatrix(
	    int i1, int i2) const
	{
	  TMVAssert(OKSubTriMatrix(i1,i2,1));
	  return GenLowerTriMatrix<T>::SubTriMatrix(i1-1,i2);
	}

	inline ConstLowerTriMatrixView<T,FortranStyle> SubTriMatrix(
	    int i1, int i2, int istep) const
	{
	  TMVAssert(OKSubTriMatrix(i1,i2,istep));
	  return GenLowerTriMatrix<T>::SubTriMatrix(i1-1,i2-1+istep,istep);
	}

	inline ConstLowerTriMatrixView<T,FortranStyle> OffDiag(
	    size_t noff=1) const
	{ return GenLowerTriMatrix<T>::OffDiag(noff); }

	inline ConstLowerTriMatrixView<T,FortranStyle> MakeUnitDiag() const
	{ return GenLowerTriMatrix<T>::MakeUnitDiag(); }

	inline ConstLowerTriMatrixView<T,FortranStyle> View() const
	{ return GenLowerTriMatrix<T>::View(); }

	inline ConstUpperTriMatrixView<T,FortranStyle> Transpose() const
	{ return GenLowerTriMatrix<T>::Transpose(); }

	inline ConstLowerTriMatrixView<T,FortranStyle> Conjugate() const
	{ return GenLowerTriMatrix<T>::Conjugate(); }

	inline ConstUpperTriMatrixView<T,FortranStyle> Adjoint() const
	{ return GenLowerTriMatrix<T>::Adjoint(); }

	inline ConstLowerTriMatrixView<RealType(T),FortranStyle> Real() const
	{ return GenLowerTriMatrix<T>::Real(); }

	inline ConstLowerTriMatrixView<RealType(T),FortranStyle> Imag() const
	{ return GenLowerTriMatrix<T>::Imag(); }

	using ConstLowerTriMatrixView<T,CStyle>::size;
	using GenLowerTriMatrix<T>::isunit;

      protected :

	using GenLowerTriMatrix<T>::okij;
	using ConstLowerTriMatrixView<T,CStyle>::cref;

      private :

	inline void operator=(const ConstLowerTriMatrixView<T,FortranStyle>&) 
	{ TMVAssert(FALSE); }

    }; // FortranStyle ConstLowerTriMatrixView

  template <class T, IndexStyle I> class UpperTriMatrixView : 
    public GenUpperTriMatrix<T>
  {

    public:

      //
      // Constructors
      //

      inline UpperTriMatrixView(const UpperTriMatrixView<T,I>& rhs) : 
	itsm(rhs.itsm), itss(rhs.itss), itssi(rhs.itssi), itssj(rhs.itssj),
	itsdiag(rhs.itsdiag), itsstor(rhs.itsstor), itsct(rhs.itsct)
	  DEFFIRSTLAST(rhs.first,rhs.last) {}

      inline UpperTriMatrixView(
	  T* _m, size_t _s, int _si, int _sj,
	  DiagType _dt, StorageType _stor, ConjType _ct
	  PARAMFIRSTLAST(T) ) :
	itsm(_m), itss(_s), itssi(_si), itssj(_sj),
	itsdiag(_dt), itsstor(_stor), itsct(_ct) DEFFIRSTLAST(_first,_last)
      {
	TMVAssert(_stor==RowMajor ? _sj == 1 : _stor==ColMajor ?
	    _si==1 : true);
      }

      virtual inline ~UpperTriMatrixView() {} 

      //
      // Op=
      //

      inline const UpperTriMatrixView<T,I>& operator=(
	  const UpperTriMatrixView<T,I>& m2) const
      {
	TMVAssert(size() == m2.size());
	m2.AssignToU(*this);
	return *this; 
      }

      inline const UpperTriMatrixView<T,I>& operator=(
	  const GenUpperTriMatrix<RealType(T)>& m2) const
      { 
	TMVAssert(size() == m2.size());
	m2.AssignToU(*this);
	return *this; 
      }

      inline const UpperTriMatrixView<T,I>& operator=(
	  const GenUpperTriMatrix<ComplexType(T)>& m2) const
      { 
	TMVAssert(size() == m2.size());
	TMVAssert(IsComplex(T()));
	m2.AssignToU(*this);
	return *this; 
      }

      template <class T2> inline const UpperTriMatrixView<T,I>& operator=(
	  const GenUpperTriMatrix<T2>& m2) const
      { 
	TMVAssert(size() == m2.size());
	TMVAssert(IsReal(T2()) || IsComplex(T()));
	TMVAssert(!isunit() || m2.isunit());
	Copy(m2,*this); 
	return *this; 
      }

      inline const UpperTriMatrixView<T,I>& operator=(T x) const 
      { TMVAssert(!isunit() || x==T(1)); return SetToIdentity(x); }

      inline const UpperTriMatrixView<T,I>& operator=(
	  const AssignableToUpperTriMatrix<RealType(T)>& m2) const
      { 
	TMVAssert(size() == m2.size());
	TMVAssert(!(m2.dt()==NonUnitDiag && dt()==UnitDiag));
	m2.AssignToU(View());
	return *this;
      }

      inline const UpperTriMatrixView<T,I>& operator=(
	  const AssignableToUpperTriMatrix<ComplexType(T)>& m2) const
      { 
	TMVAssert(size() == m2.size());
	TMVAssert(IsComplex(T()));
	TMVAssert(!(m2.dt()==NonUnitDiag && dt()==UnitDiag));
	m2.AssignToU(View());
	return *this;
      }

      //
      // Access
      //

      inline RefType(T) operator()(size_t i,size_t j) const 
      {
	TMVAssert(i<size());
	TMVAssert(j<size());
	TMVAssert(okij(i,j));
	return ref(i,j); 
      }

      inline VectorView<T,I> row(size_t i, size_t j1, size_t j2) const 
      { 
	TMVAssert(i<size());
	TMVAssert(j1<=j2);
	TMVAssert(j2<=size());
	TMVAssert(j1==j2 || okij(i,j1));
	return VectorView<T,I>(
	    ptr()+i*stepi()+j1*stepj(),j2-j1,stepj(),ct() FIRSTLAST); 
      }

      inline VectorView<T,I> col(size_t j, size_t i1, size_t i2) const
      {
	TMVAssert(j<size());
	TMVAssert(i1<=i2);
	TMVAssert(i2<=size());
	TMVAssert(i1==i2 || okij(i2-1,j));
	return VectorView<T,I>(
	    ptr()+i1*stepi()+j*stepj(),i2-i1,stepi(),ct() FIRSTLAST); 
      }

      inline VectorView<T,I> diag() const
      {
	TMVAssert(!isunit());
	return VectorView<T,I>(ptr(),size(),stepi()+stepj(),ct() FIRSTLAST); 
      }

      inline VectorView<T,I> diag(int i) const
      {
	TMVAssert(i<=int(size())); 
	TMVAssert(isunit() ? i>0 : i>=0);
	return VectorView<T,I>(ptr()+i*stepj(),size()-i,(stepi()+stepj()),ct() 
	    FIRSTLAST);
      }

      inline VectorView<T,I> diag(int i, size_t j1, size_t j2) const
      {
	TMVAssert(isunit() ? i>0 : i>=0);
	TMVAssert(i<=int(size())); 
	TMVAssert(j1 <= j2);
	TMVAssert(j2 <= size()-i);
	const int ds = stepi()+stepj();
	return VectorView<T,I>(ptr()+i*stepj()+j1*ds,j2-j1,ds,ct() FIRSTLAST);
      }

      //
      // Modifying Functions
      //

      inline const UpperTriMatrixView<T,I>& Zero() const 
      { return SetAllTo(T(0)); }

      const UpperTriMatrixView<T,I>& SetAllTo(T x) const;

      const UpperTriMatrixView<T,I>& Clip(RealType(T) thresh) const;

      const UpperTriMatrixView<T,I>& ConjugateSelf() const;

      const UpperTriMatrixView<T,I>& InvertSelf() const;

      const UpperTriMatrixView<T,I>& SetToIdentity(T x=T(1)) const;

      //
      // SubMatrix
      //

      using GenUpperTriMatrix<T>::OKSubMatrix;
      using GenUpperTriMatrix<T>::OKSubVector;
      using GenUpperTriMatrix<T>::OKSubTriMatrix;

      inline MatrixView<T,I> SubMatrix(int i1, int i2, int j1, int j2) const
      {
	TMVAssert(OKSubMatrix(i1,i2,j1,j2,1,1));
	return MatrixView<T,I>(ptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, stepi(),stepj(),stor(),ct() FIRSTLAST );
      }

      inline MatrixView<T,I> SubMatrix(int i1, int i2, int j1, int j2,
	  int istep, int jstep) const
      {
	const StorageType newstor =
	  iscm() ? (istep == 1 ? ColMajor : NoMajor) :
	  isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
	TMVAssert(OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	return MatrixView<T,I>(ptr()+i1*stepi()+j1*stepj(),
	    (i2-i1)/istep, (j2-j1)/jstep, istep*stepi(), jstep*stepj(),
	    newstor,ct() FIRSTLAST );
      }

      inline VectorView<T,I> SubVector(int i, int j,
	  int istep, int jstep, size_t size) const
      {
	TMVAssert(OKSubVector(i,j,istep,jstep,size));
	return VectorView<T,I>(ptr()+i*stepi()+j*stepj(),size,
	    istep*stepi()+jstep*stepj(),ct() FIRSTLAST );
      }

      inline UpperTriMatrixView<T,I> SubTriMatrix(int i1, int i2) const
      {
	TMVAssert(OKSubTriMatrix(i1,i2,1));
	return UpperTriMatrixView<T,I>(ptr()+i1*(stepi()+stepj()),i2-i1,
	    stepi(),stepj(),dt(),stor(),ct() FIRSTLAST);
      }

      inline UpperTriMatrixView<T,I> SubTriMatrix(int i1, int i2,
	  int istep) const
      {
	TMVAssert(OKSubTriMatrix(i1,i2,istep));
	return UpperTriMatrixView<T,I>(ptr()+i1*(stepi()+stepj()),
	    (i2-i1)/istep,istep*stepi(),istep*stepj(),dt(),
	    istep==1 ? stor() : NoMajor,ct() FIRSTLAST);
      }

      inline UpperTriMatrixView<T,I> OffDiag(size_t noff=1) const
      {
	TMVAssert(size() >= noff);
	return UpperTriMatrixView<T,I>(ptr()+noff*stepj(),size()-noff,
	    stepi(),stepj(),NonUnitDiag,stor(),ct() FIRSTLAST);
      }

      inline UpperTriMatrixView<T,I> MakeUnitDiag() const
      { 
	return UpperTriMatrixView<T,I>(ptr(),size(),
	    stepi(),stepj(),UnitDiag,stor(),ct() FIRSTLAST);
      }

      inline UpperTriMatrixView<T,I> View() const
      { return *this; }

      inline LowerTriMatrixView<T,I> Transpose() const
      {
	return LowerTriMatrixView<T,I>(ptr(),size(),
	    stepj(),stepi(),dt(),TransOf(stor()),ct() FIRSTLAST);
      }

      inline UpperTriMatrixView<T,I> Conjugate() const
      {
	return UpperTriMatrixView<T,I>(ptr(),size(),
	    stepi(),stepj(),dt(),stor(),ConjOf(T,ct()) FIRSTLAST);
      }

      inline LowerTriMatrixView<T,I> Adjoint() const
      {
	return LowerTriMatrixView<T,I>(ptr(),size(),
	    stepj(),stepi(),dt(),TransOf(stor()),ConjOf(T,ct()) 
	    FIRSTLAST);
      }

      inline UpperTriMatrixView<RealType(T)> Real() const
      {
	return UpperTriMatrixView<RealType(T)>(
	    reinterpret_cast<RealType(T)*>(ptr()),size(),
	    IsReal(T()) ? stepi() : 2*stepi(),
	    IsReal(T()) ? stepj() : 2*stepj(),
	    dt(), IsReal(T()) ? stor() : NoMajor, NonConj
#ifdef TMVFLDEBUG
	    ,reinterpret_cast<const RealType(T)*>(first)
	    ,reinterpret_cast<const RealType(T)*>(last)
#endif
	    );
      }

      inline UpperTriMatrixView<RealType(T)> Imag() const
      {
	TMVAssert(IsComplex(T()));
	TMVAssert(!isunit());
	return UpperTriMatrixView<RealType(T)>(
	    reinterpret_cast<RealType(T)*>(ptr())+1,
	    size(),2*stepi(),2*stepj(),dt(),NoMajor, NonConj
#ifdef TMVFLDEBUG
	    ,reinterpret_cast<const RealType(T)*>(first)+1
	    ,reinterpret_cast<const RealType(T)*>(last)+1
#endif
	    );
      }

      //
      // I/O
      //

      void Read(std::istream& is) const;

      inline size_t size() const { return itss; }
      inline const T* cptr() const { return itsm; }
      inline T* ptr() const { return itsm; }
      inline int stepi() const { return itssi; }
      inline int stepj() const { return itssj; }
      using GenUpperTriMatrix<T>::isconj;
      using GenUpperTriMatrix<T>::isrm;
      using GenUpperTriMatrix<T>::iscm;
      using GenUpperTriMatrix<T>::isunit;
      inline StorageType stor() const { return itsstor; }
      inline DiagType dt() const { return itsdiag; }
      inline ConjType ct() const { return itsct; }

#if 0
      // This makes it easier for some things to compile
      // The ones that work are overridden
      template <class T2> inline void operator=(
	  const BaseMatrix<T2>&) const { TMVAssert(FALSE); }
      template <class T2> inline void operator+=(
	  const BaseMatrix<T2>&) const { TMVAssert(FALSE); }
      template <class T2> inline void operator-=(
	  const BaseMatrix<T2>&) const { TMVAssert(FALSE); }
      template <class T2> inline void operator*=(
	  const BaseMatrix<T2>&) const { TMVAssert(FALSE); }
      template <class T2> inline void operator/=(
	  const BaseMatrix<T2>&) const { TMVAssert(FALSE); }
      template <class T2> inline void operator%=(
	  const BaseMatrix<T2>&) const { TMVAssert(FALSE); }
#endif

    protected :

      T*const itsm;
      const size_t itss;
      const int itssi;
      const int itssj;
      DiagType itsdiag;
      StorageType itsstor;
      ConjType itsct;

#ifdef TMVFLDEBUG
    public:
      const T*const first;
      const T*const last;
    protected:
#endif

      using GenUpperTriMatrix<T>::okij;
      RefType(T) ref(size_t i, size_t j) const;

  }; // UpperTriMatrixView

  template <class T, IndexStyle I> class LowerTriMatrixView : 
    public GenLowerTriMatrix<T>
  {

    public:

      //
      // Constructors
      //

      inline LowerTriMatrixView(const LowerTriMatrixView<T,I>& rhs) : 
	itsm(rhs.itsm), itss(rhs.itss), itssi(rhs.itssi), itssj(rhs.itssj),
	itsdiag(rhs.itsdiag), itsstor(rhs.itsstor), itsct(rhs.itsct)
	  DEFFIRSTLAST(rhs.first,rhs.last) {}

      inline LowerTriMatrixView(
	  T* _m, size_t _s, int _si, int _sj,
	  DiagType _dt, StorageType _stor, ConjType _ct
	  PARAMFIRSTLAST(T) ) :
	itsm(_m), itss(_s), itssi(_si), itssj(_sj),
	itsdiag(_dt), itsstor(_stor), itsct(_ct) DEFFIRSTLAST(_first,_last)
      {
	TMVAssert(_stor==RowMajor ? _sj == 1 : _stor==ColMajor ?
	    _si==1 : true);
      }

      virtual inline ~LowerTriMatrixView() {} 

      //
      // Op=
      //

      inline const LowerTriMatrixView<T,I>& operator=(
	  const LowerTriMatrixView<T,I>& m2) const
      {
	TMVAssert(size() == m2.size());
	m2.AssignToL(*this);
	return *this; 
      }

      inline const LowerTriMatrixView<T,I>& operator=(
	  const GenLowerTriMatrix<RealType(T)>& m2) const
      {
	TMVAssert(size() == m2.size());
	m2.AssignToL(*this);
	return *this; 
      }

      inline const LowerTriMatrixView<T,I>& operator=(
	  const GenLowerTriMatrix<ComplexType(T)>& m2) const
      {
	TMVAssert(size() == m2.size());
	TMVAssert(IsComplex(T()));
	m2.AssignToL(*this);
	return *this; 
      }

      template <class T2> inline const LowerTriMatrixView<T,I>& operator=(
	  const GenLowerTriMatrix<T2>& m2) const
      { 
	TMVAssert(size() == m2.size());
	TMVAssert(!isunit() || m2.isunit());
	Copy(m2.Transpose(),Transpose()); 
	return *this; 
      }

      inline const LowerTriMatrixView<T,I>& operator=(T x) const 
      { 
	TMVAssert(!isunit()); 
	return SetToIdentity(x); 
      }

      inline const LowerTriMatrixView<T,I>& operator=(
	  const AssignableToLowerTriMatrix<RealType(T)>& m2) const
      { 
	TMVAssert(size() == m2.size());
	TMVAssert(!(m2.dt()==NonUnitDiag && dt()==UnitDiag));
	m2.AssignToL(View());
	return *this;
      }

      inline const LowerTriMatrixView<T,I>& operator=(
	  const AssignableToLowerTriMatrix<ComplexType(T)>& m2) const
      { 
	TMVAssert(size() == m2.size());
	TMVAssert(IsComplex(T()));
	TMVAssert(!(m2.dt()==NonUnitDiag && dt()==UnitDiag));
	m2.AssignToL(View());
	return *this;
      }

      //
      // Access
      //

      inline RefType(T) operator()(size_t i,size_t j) const 
      {
	TMVAssert(i<size());
	TMVAssert(j<size());
	TMVAssert(okij(i,j));
	return ref(i,j); 
      }

      inline VectorView<T,I> row(size_t i, size_t j1, size_t j2) const 
      { 
	TMVAssert(i<size());
	TMVAssert(j1<=j2);
	TMVAssert(j2<=size());
	TMVAssert(j1==j2 || okij(i,j2-1));
	return VectorView<T,I>(ptr()+i*stepi()+j1*stepj(),j2-j1,stepj(),ct() 
	    FIRSTLAST); 
      }

      inline VectorView<T,I> col(size_t j, size_t i1, size_t i2) const
      {
	TMVAssert(j<size());
	TMVAssert(i1<=i2);
	TMVAssert(i2<=size());
	TMVAssert(i1==i2 || okij(i1,j));
	return VectorView<T,I>(ptr()+i1*stepi()+j*stepj(),i2-i1,stepi(),ct() 
	    FIRSTLAST); 
      }

      inline VectorView<T,I> diag() const
      {
	TMVAssert(!isunit());
	return VectorView<T,I>(ptr(),size(),stepi()+stepj(),ct() FIRSTLAST); 
      }

      inline VectorView<T,I> diag(int i) const
      {
	TMVAssert(-i<=int(size())); 
	TMVAssert(isunit() ? i<0 : i<=0);
	return VectorView<T,I>(ptr()-i*stepi(),size()+i,stepi()+stepj(),ct() 
	    FIRSTLAST);
      }

      inline VectorView<T,I> diag(int i, size_t j1, size_t j2) const
      {
	TMVAssert(isunit() ? i<0 : i<=0);
	TMVAssert(-i<=int(size())); 
	TMVAssert(j1 <= j2);
	TMVAssert(j2 <= size()+i);
	const int ds = stepi()+stepj();
	return VectorView<T,I>(ptr()-i*stepi()+j1*ds,j2-j1,ds,ct() FIRSTLAST);
      }

      //
      // Modifying Functions
      //

      inline const LowerTriMatrixView<T,I>& Zero() const 
      { return SetAllTo(T(0)); }

      inline const LowerTriMatrixView<T,I>& SetAllTo(T x) const
      { Transpose().SetAllTo(x); return *this; }

      inline const LowerTriMatrixView<T,I>& Clip(RealType(T) thresh) const
      { Transpose().Clip(thresh); return *this; }

      inline const LowerTriMatrixView<T,I>& ConjugateSelf() const
      { Transpose().ConjugateSelf(); return *this; }

      inline const LowerTriMatrixView<T,I>& InvertSelf() const
      { Transpose().InvertSelf(); return *this; }

      inline const LowerTriMatrixView<T,I>& SetToIdentity(T x=T(1)) const
      { Transpose().SetToIdentity(x); return *this; }

      //
      // SubMatrix
      //

      using GenLowerTriMatrix<T>::OKSubMatrix;
      using GenLowerTriMatrix<T>::OKSubVector;
      using GenLowerTriMatrix<T>::OKSubTriMatrix;

      inline MatrixView<T,I> SubMatrix(int i1, int i2, int j1, int j2) const
      {
	TMVAssert(OKSubMatrix(i1,i2,j1,j2,1,1));
	return MatrixView<T,I>(ptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, stepi(),stepj(),stor(),ct() FIRSTLAST );
      }

      inline MatrixView<T,I> SubMatrix(int i1, int i2, int j1, int j2,
	  int istep, int jstep) const
      {
	const StorageType newstor =
	  iscm() ? (istep == 1 ? ColMajor : NoMajor) :
	  isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
	TMVAssert(OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	return MatrixView<T,I>(ptr()+i1*stepi()+j1*stepj(),
	    (i2-i1)/istep, (j2-j1)/jstep, istep*stepi(), jstep*stepj(),
	    newstor,ct() FIRSTLAST );
      }

      inline VectorView<T,I> SubVector(int i, int j,
	  int istep, int jstep, size_t size) const
      {
	TMVAssert(OKSubVector(i,j,istep,jstep,size));
	return VectorView<T,I>(ptr()+i*stepi()+j*stepj(),size,
	    istep*stepi()+jstep*stepj(),ct() FIRSTLAST );
      }

      inline LowerTriMatrixView<T,I> SubTriMatrix(int i1, int i2) const
      {
	TMVAssert(OKSubTriMatrix(i1,i2,1));
	return LowerTriMatrixView<T,I>(ptr()+i1*(stepi()+stepj()),
	    i2-i1,stepi(),stepj(),dt(),stor(),ct() FIRSTLAST);
      }

      inline LowerTriMatrixView<T,I> SubTriMatrix(int i1, int i2,
	  int istep) const
      {
	TMVAssert(OKSubTriMatrix(i1,i2,istep));
	return LowerTriMatrixView<T,I>(ptr()+i1*(stepi()+stepj()),
	    (i2-i1)/istep,istep*stepi(),istep*stepj(),dt(),
	    istep==1 ? stor() : NoMajor,ct() FIRSTLAST);
      }

      inline LowerTriMatrixView<T,I> OffDiag(size_t noff=1) const
      {
	TMVAssert(size() >= noff);
	return LowerTriMatrixView<T,I>(ptr()+noff*stepi(),size()-noff,
	    stepi(),stepj(),NonUnitDiag,stor(),ct() FIRSTLAST);
      }

      inline LowerTriMatrixView<T,I> MakeUnitDiag() const
      { 
	return LowerTriMatrixView<T,I>(ptr(),size(),
	    stepi(),stepj(),UnitDiag,stor(),ct() FIRSTLAST);
      }

      inline LowerTriMatrixView<T,I> View() const
      { return *this; }

      inline UpperTriMatrixView<T,I> Transpose() const
      {
	return UpperTriMatrixView<T,I>(ptr(),size(),
	    stepj(),stepi(),dt(),TransOf(stor()),ct() FIRSTLAST);
      }

      inline LowerTriMatrixView<T,I> Conjugate() const
      {
	return LowerTriMatrixView<T,I>(ptr(),size(),
	    stepi(),stepj(),dt(),stor(),ConjOf(T,ct()) FIRSTLAST);
      }

      inline UpperTriMatrixView<T,I> Adjoint() const
      {
	return UpperTriMatrixView<T,I>(ptr(),size(),
	    stepj(),stepi(),dt(),TransOf(stor()),ConjOf(T,ct()) 
	    FIRSTLAST);
      }

      inline LowerTriMatrixView<RealType(T)> Real() const
      {
	return LowerTriMatrixView<RealType(T)>(
	    reinterpret_cast<RealType(T)*>(ptr()),size(),
	    IsReal(T()) ? stepi() : 2*stepi(),
	    IsReal(T()) ? stepj() : 2*stepj(),
	    dt(), IsReal(T()) ? stor() : NoMajor,NonConj
#ifdef TMVFLDEBUG
	    ,reinterpret_cast<const RealType(T)*>(first)
	    ,reinterpret_cast<const RealType(T)*>(last)
#endif
	    );
      }

      inline LowerTriMatrixView<RealType(T)> Imag() const
      {
	TMVAssert(IsComplex(T()));
	TMVAssert(!isunit());
	return LowerTriMatrixView<RealType(T)>(
	    reinterpret_cast<RealType(T)*>(ptr())+1,
	    size(),2*stepi(),2*stepj(),dt(),NoMajor,NonConj
#ifdef TMVFLDEBUG
	    ,reinterpret_cast<const RealType(T)*>(first)+1
	    ,reinterpret_cast<const RealType(T)*>(last)+1
#endif
	    );
      }

      //
      // I/O
      //

      void Read(std::istream& is) const;

      inline size_t size() const { return itss; }
      inline const T* cptr() const { return itsm; }
      inline T* ptr() const { return itsm; }
      inline int stepi() const { return itssi; }
      inline int stepj() const { return itssj; }
      using GenLowerTriMatrix<T>::isconj;
      using GenLowerTriMatrix<T>::isrm;
      using GenLowerTriMatrix<T>::iscm;
      using GenLowerTriMatrix<T>::isunit;
      inline StorageType stor() const { return itsstor; }
      inline DiagType dt() const { return itsdiag; }
      inline ConjType ct() const { return itsct; }

#if 0
      // This makes it easier for some things to compile
      // The ones that work are overridden
      template <class T2> inline void operator=(
	  const BaseMatrix<T2>&) const { TMVAssert(FALSE); }
      template <class T2> inline void operator+=(
	  const BaseMatrix<T2>&) const { TMVAssert(FALSE); }
      template <class T2> inline void operator-=(
	  const BaseMatrix<T2>&) const { TMVAssert(FALSE); }
      template <class T2> inline void operator*=(
	  const BaseMatrix<T2>&) const { TMVAssert(FALSE); }
      template <class T2> inline void operator/=(
	  const BaseMatrix<T2>&) const { TMVAssert(FALSE); }
      template <class T2> inline void operator%=(
	  const BaseMatrix<T2>&) const { TMVAssert(FALSE); }
#endif

    protected :

      T*const itsm;
      const size_t itss;
      const int itssi;
      const int itssj;
      DiagType itsdiag;
      StorageType itsstor;
      ConjType itsct;

#ifdef TMVFLDEBUG
    public:
      const T*const first;
      const T*const last;
    protected:
#endif

      using GenLowerTriMatrix<T>::okij;
      RefType(T) ref(size_t i, size_t j) const;

  }; // LowerTriMatrixView

  template <class T> class UpperTriMatrixView<T,FortranStyle> : 
    public UpperTriMatrixView<T,CStyle>
    {

      public:

	//
	// Constructors
	//

	inline UpperTriMatrixView(
	    const UpperTriMatrixView<T,FortranStyle>& rhs) : 
	  UpperTriMatrixView<T,CStyle>(rhs) {}

	inline UpperTriMatrixView(const UpperTriMatrixView<T,CStyle>& rhs) : 
	  UpperTriMatrixView<T,CStyle>(rhs) {}

	inline UpperTriMatrixView(
	    T* _m, size_t _s, int _si, int _sj,
	    DiagType indt, StorageType instor, ConjType inct 
	    PARAMFIRSTLAST(T) ) :
	  UpperTriMatrixView<T,CStyle>(_m,_s,_si,_sj,indt,instor,inct 
	      FIRSTLAST1(_first,_last) ) {}

	virtual inline ~UpperTriMatrixView() {} 

	//
	// Op=
	//

	inline const UpperTriMatrixView<T,FortranStyle>& operator=(
	    const UpperTriMatrixView<T,FortranStyle>& m2) const
	{ UpperTriMatrixView<T,CStyle>::operator=(m2); return *this; }

	inline const UpperTriMatrixView<T,FortranStyle>& operator=(
	    const GenUpperTriMatrix<RealType(T)>& m2) const
	{ UpperTriMatrixView<T,CStyle>::operator=(m2); return *this; }

	inline const UpperTriMatrixView<T,FortranStyle>& operator=(
	    const GenUpperTriMatrix<ComplexType(T)>& m2) const
	{ UpperTriMatrixView<T,CStyle>::operator=(m2); return *this; }

	template <class T2> 
	  inline const UpperTriMatrixView<T,FortranStyle>& operator=(
	      const GenUpperTriMatrix<T2>& m2) const
	  { UpperTriMatrixView<T,CStyle>::operator=(m2); return *this; }

	inline const UpperTriMatrixView<T,FortranStyle>& operator=(T x) const 
	{ UpperTriMatrixView<T,CStyle>::operator=(x); return *this; }

	inline const UpperTriMatrixView<T,FortranStyle>& operator=(
	    const AssignableToUpperTriMatrix<RealType(T)>& m2) const
	{ UpperTriMatrixView<T,CStyle>::operator=(m2); return *this; }

	inline const UpperTriMatrixView<T,FortranStyle>& operator=(
	    const AssignableToUpperTriMatrix<ComplexType(T)>& m2) const
	{ UpperTriMatrixView<T,CStyle>::operator=(m2); return *this; }

	//
	// Access
	//

	inline RefType(T) operator()(size_t i,size_t j) const 
	{ 
	  TMVAssert(i>0 && i <= size());
	  TMVAssert(j>0 && j <= size());
	  TMVAssert(okij(i-1,j-1));
	  return ref(i-1,j-1); 
	}

	inline VectorView<T,FortranStyle> row(
	    size_t i, size_t j1, size_t j2) const 
	{ 
	  TMVAssert(i>0 && i<=size());
	  TMVAssert(j1>0 && j1<=j2 && j2<=size());
	  return UpperTriMatrixView<T,CStyle>::row(i-1,j1-1,j2);
	}

	inline VectorView<T,FortranStyle> col(
	    size_t j, size_t i1, size_t i2) const
	{
	  TMVAssert(j>0 && j<=size());
	  TMVAssert(i1>0 && i1<=i2 && i2<=size());
	  return UpperTriMatrixView<T,CStyle>::col(j-1,i1-1,i2);
	}

	inline VectorView<T,FortranStyle> diag() const
	{ return UpperTriMatrixView<T,CStyle>::diag(); }

	inline VectorView<T,FortranStyle> diag(int i) const
	{ return UpperTriMatrixView<T,CStyle>::diag(i); }

	inline VectorView<T,FortranStyle> diag(
	    int i, size_t j1, size_t j2) const
	{
	  TMVAssert(j1>0);
	  return UpperTriMatrixView<T,CStyle>::diag(i,j1-1,j2); 
	}

	//
	// Modifying Functions
	//

	inline const UpperTriMatrixView<T,FortranStyle>& Zero() const 
	{ UpperTriMatrixView<T,CStyle>::Zero(); return *this; }

	inline const UpperTriMatrixView<T,FortranStyle>& SetAllTo(T x) const
	{ UpperTriMatrixView<T,CStyle>::SetAllTo(x); return *this; }

	inline const UpperTriMatrixView<T,FortranStyle>& Clip(
	    RealType(T) thresh) const
	{ UpperTriMatrixView<T,CStyle>::Clip(thresh); return *this; }

	inline const UpperTriMatrixView<T,FortranStyle>& ConjugateSelf() const
	{ UpperTriMatrixView<T,CStyle>::ConjugateSelf(); return *this; }

	inline const LowerTriMatrixView<T,FortranStyle>& InvertSelf() const
	{ UpperTriMatrixView<T,CStyle>::InvertSelf(); return *this; }

	inline const UpperTriMatrixView<T,FortranStyle>& SetToIdentity(
	    T x=T(1)) const
	{ UpperTriMatrixView<T,CStyle>::SetToIdentity(x); return *this; }

	//
	// SubMatrix
	//

	inline bool OKSubMatrix(
	    int i1, int i2, int j1, int j2, int istep, int jstep) const
	{
	  return ConstUpperTriMatrixView<T,FortranStyle>(*this).OKSubMatrix(
	      i1,i2,j1,j2,istep,jstep); 
	}

	inline bool OKSubVector(
	    int i, int j, int istep, int jstep, size_t s) const
	{
	  return ConstUpperTriMatrixView<T,FortranStyle>(*this).OKSubVector(
	      i,j,istep,jstep,s);
	}

	inline bool OKSubTriMatrix(int i1, int i2, int istep) const
	{
	  return ConstUpperTriMatrixView<T,FortranStyle>(*this).OKSubTriMatrix(
	      i1,i2,istep);
	}

	inline MatrixView<T,FortranStyle> SubMatrix(
	    int i1, int i2, int j1, int j2) const
	{
	  TMVAssert(OKSubMatrix(i1,i2,j1,j2,1,1));
	  return UpperTriMatrixView<T,CStyle>::SubMatrix(i1-1,i2,j1-1,j2);
	}

	inline MatrixView<T,FortranStyle> SubMatrix(
	    int i1, int i2, int j1, int j2, int istep, int jstep) const
	{
	  TMVAssert(OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	  return UpperTriMatrixView<T,CStyle>::SubMatrix(
	      i1-1,i2-1+istep,j1-1,j2-1+jstep,istep,jstep);
	}

	inline VectorView<T,FortranStyle> SubVector(size_t i, size_t j,
	    int istep, int jstep, size_t s) const
	{
	  TMVAssert(OKSubVector(i,j,istep,jstep,s));
	  return UpperTriMatrixView<T,CStyle>::SubVector(
	      i-1,j-1,istep,jstep,s);
	}

	inline UpperTriMatrixView<T,FortranStyle> SubTriMatrix(
	    int i1, int i2) const
	{
	  TMVAssert(OKSubTriMatrix(i1,i2,1));
	  return UpperTriMatrixView<T,CStyle>::SubTriMatrix(i1-1,i2);
	}

	inline UpperTriMatrixView<T,FortranStyle> SubTriMatrix(int i1, int i2,
	    int istep) const
	{
	  TMVAssert(OKSubTriMatrix(i1,i2,istep));
	  return UpperTriMatrixView<T,CStyle>::SubTriMatrix(
	      i1-1,i2-1+istep,istep);
	}

	inline UpperTriMatrixView<T,FortranStyle> OffDiag(size_t noff=1) const
	{
	  TMVAssert(size() >= noff);
	  return UpperTriMatrixView<T,CStyle>::OffDiag(noff);
	}

	inline UpperTriMatrixView<T,FortranStyle> MakeUnitDiag() const
	{ return UpperTriMatrixView<T,CStyle>::MakeUnitDiag(); }

	inline UpperTriMatrixView<T,FortranStyle> View() const
	{ return *this; }

	inline LowerTriMatrixView<T,FortranStyle> Transpose() const
	{ return UpperTriMatrixView<T,CStyle>::Transpose(); }

	inline UpperTriMatrixView<T,FortranStyle> Conjugate() const
	{ return UpperTriMatrixView<T,CStyle>::Conjugate(); }

	inline LowerTriMatrixView<T,FortranStyle> Adjoint() const
	{ return UpperTriMatrixView<T,CStyle>::Adjoint(); }

	inline UpperTriMatrixView<RealType(T)> Real() const
	{ return UpperTriMatrixView<T,CStyle>::Real(); }

	inline UpperTriMatrixView<RealType(T)> Imag() const
	{ return UpperTriMatrixView<T,CStyle>::Imag(); }

	template <class T2> inline void operator=(const BaseMatrix<T2>&) const 
	{ TMVAssert(FALSE); }

	using UpperTriMatrixView<T,CStyle>::size;

      protected :

	using GenUpperTriMatrix<T>::okij;
	using UpperTriMatrixView<T,CStyle>::ref;

    }; // FortranStyle UpperTriMatrixView

  template <class T> class LowerTriMatrixView<T,FortranStyle> : 
    public LowerTriMatrixView<T,CStyle>
    {

      public:

	//
	// Constructors
	//

	inline LowerTriMatrixView(
	    const LowerTriMatrixView<T,FortranStyle>& rhs) : 
	  LowerTriMatrixView<T,CStyle>(rhs) {}

	inline LowerTriMatrixView(const LowerTriMatrixView<T,CStyle>& rhs) : 
	  LowerTriMatrixView<T,CStyle>(rhs) {}

	inline LowerTriMatrixView(
	    T* _m, size_t _s, int _si, int _sj,
	    DiagType indt, StorageType instor, ConjType inct 
	    PARAMFIRSTLAST(T) ) :
	  LowerTriMatrixView<T,CStyle>(_m,_s,_si,_sj,indt,instor,inct 
	      FIRSTLAST1(_first,_last) ) {}

	virtual inline ~LowerTriMatrixView() {} 

	//
	// Op=
	//

	inline const LowerTriMatrixView<T,FortranStyle>& operator=(
	    const LowerTriMatrixView<T,FortranStyle>& m2) const
	{ LowerTriMatrixView<T,CStyle>::operator=(m2); return *this; }

	inline const LowerTriMatrixView<T,FortranStyle>& operator=(
	    const GenLowerTriMatrix<RealType(T)>& m2) const
	{ LowerTriMatrixView<T,CStyle>::operator=(m2); return *this; }

	inline const LowerTriMatrixView<T,FortranStyle>& operator=(
	    const GenLowerTriMatrix<ComplexType(T)>& m2) const
	{ LowerTriMatrixView<T,CStyle>::operator=(m2); return *this; }

	template <class T2> 
	  inline const LowerTriMatrixView<T,FortranStyle>& operator=(
	      const GenLowerTriMatrix<T2>& m2) const
	  { LowerTriMatrixView<T,CStyle>::operator=(m2); return *this; }

	inline const LowerTriMatrixView<T,FortranStyle>& operator=(T x) const 
	{ LowerTriMatrixView<T,CStyle>::operator=(x); return *this; }

	inline const LowerTriMatrixView<T,FortranStyle>& operator=(
	    const AssignableToLowerTriMatrix<RealType(T)>& m2) const
	{ LowerTriMatrixView<T,CStyle>::operator=(m2); return *this; }

	inline const LowerTriMatrixView<T,FortranStyle>& operator=(
	    const AssignableToLowerTriMatrix<ComplexType(T)>& m2) const
	{ LowerTriMatrixView<T,CStyle>::operator=(m2); return *this; }

	//
	// Access
	//

	inline RefType(T) operator()(size_t i,size_t j) const 
	{ 
	  TMVAssert(i>0 && i <= size());
	  TMVAssert(j>0 && j <= size());
	  TMVAssert(okij(i-1,j-1));
	  return ref(i-1,j-1); 
	}

	inline VectorView<T,FortranStyle> row(
	    size_t i, size_t j1, size_t j2) const 
	{ 
	  TMVAssert(i>0 && i<=size());
	  TMVAssert(j1>0 && j1<=j2 && j2<=size());
	  return LowerTriMatrixView<T,CStyle>::row(i-1,j1-1,j2);
	}

	inline VectorView<T,FortranStyle> col(
	    size_t j, size_t i1, size_t i2) const
	{
	  TMVAssert(j>0 && j<=size());
	  TMVAssert(i1>0 && i1<=i2 && i2<=size());
	  return LowerTriMatrixView<T,CStyle>::col(j-1,i1-1,i2);
	}

	inline VectorView<T,FortranStyle> diag() const
	{ return LowerTriMatrixView<T,CStyle>::diag(); }

	inline VectorView<T,FortranStyle> diag(int i) const
	{ return LowerTriMatrixView<T,CStyle>::diag(i); }

	inline VectorView<T,FortranStyle> diag(
	    int i, size_t j1, size_t j2) const
	{
	  TMVAssert(j1>0);
	  return LowerTriMatrixView<T,CStyle>::diag(i,j1-1,j2); 
	}

	//
	// Modifying Functions
	//

	inline const LowerTriMatrixView<T,FortranStyle>& Zero() const 
	{ LowerTriMatrixView<T,CStyle>::Zero(); return *this; }

	inline const LowerTriMatrixView<T,FortranStyle>& SetAllTo(T x) const
	{ LowerTriMatrixView<T,CStyle>::SetAllTo(x); return *this; }

	inline const LowerTriMatrixView<T,FortranStyle>& Clip(
	    RealType(T) thresh) const
	{ LowerTriMatrixView<T,CStyle>::Clip(thresh); return *this; }

	inline const LowerTriMatrixView<T,FortranStyle>& ConjugateSelf() const
	{ LowerTriMatrixView<T,CStyle>::ConjugateSelf(); return *this; }

	inline const LowerTriMatrixView<T,FortranStyle>& InvertSelf() const
	{ LowerTriMatrixView<T,CStyle>::InvertSelf(); return *this; }

	inline const LowerTriMatrixView<T,FortranStyle>& SetToIdentity(
	    T x=T(1)) const
	{ LowerTriMatrixView<T,CStyle>::SetToIdentity(x); return *this; }

	//
	// SubMatrix
	//

	inline bool OKSubMatrix(
	    int i1, int i2, int j1, int j2, int istep, int jstep) const
	{
	  return ConstLowerTriMatrixView<T,FortranStyle>(*this).OKSubMatrix(
	      i1,i2,j1,j2,istep,jstep); 
	}

	inline bool OKSubVector(
	    int i, int j, int istep, int jstep, size_t s) const
	{
	  return ConstLowerTriMatrixView<T,FortranStyle>(*this).OKSubVector(
	      i,j,istep,jstep,s);
	}

	inline bool OKSubTriMatrix(int i1, int i2, int istep) const
	{
	  return ConstLowerTriMatrixView<T,FortranStyle>(*this).OKSubTriMatrix(
	      i1,i2,istep);
	}

	inline MatrixView<T,FortranStyle> SubMatrix(
	    int i1, int i2, int j1, int j2) const
	{
	  TMVAssert(OKSubMatrix(i1,i2,j1,j2,1,1));
	  return LowerTriMatrixView<T,CStyle>::SubMatrix(i1-1,i2,j1-1,j2);
	}

	inline MatrixView<T,FortranStyle> SubMatrix(
	    int i1, int i2, int j1, int j2, int istep, int jstep) const
	{
	  TMVAssert(OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	  return LowerTriMatrixView<T,CStyle>::SubMatrix(
	      i1-1,i2-1+istep,j1-1,j2-1+jstep,istep,jstep);
	}

	inline VectorView<T,FortranStyle> SubVector(int i, int j,
	    int istep, int jstep, size_t s) const
	{
	  TMVAssert(OKSubVector(i,j,istep,jstep,s));
	  return LowerTriMatrixView<T,CStyle>::SubVector(
	      i-1,j-1,istep,jstep,s);
	}

	inline LowerTriMatrixView<T,FortranStyle> SubTriMatrix(
	    int i1, int i2) const
	{
	  TMVAssert(OKSubTriMatrix(i1,i2,1));
	  return LowerTriMatrixView<T,CStyle>::SubTriMatrix(i1-1,i2);
	}

	inline LowerTriMatrixView<T,FortranStyle> SubTriMatrix(int i1, int i2,
	    int istep) const
	{
	  TMVAssert(OKSubTriMatrix(i1,i2,istep));
	  return LowerTriMatrixView<T,CStyle>::SubTriMatrix(
	      i1-1,i2-1+istep,istep);
	}

	inline LowerTriMatrixView<T,FortranStyle> OffDiag(size_t noff=1) const
	{
	  TMVAssert(size() >= noff);
	  return LowerTriMatrixView<T,CStyle>::OffDiag(noff);
	}

	inline LowerTriMatrixView<T,FortranStyle> MakeUnitDiag() const
	{ return LowerTriMatrixView<T,CStyle>::MakeUnitDiag(); }

	inline LowerTriMatrixView<T,FortranStyle> View() const
	{ return *this; }

	inline UpperTriMatrixView<T,FortranStyle> Transpose() const
	{ return LowerTriMatrixView<T,CStyle>::Transpose(); }

	inline LowerTriMatrixView<T,FortranStyle> Conjugate() const
	{ return LowerTriMatrixView<T,CStyle>::Conjugate(); }

	inline UpperTriMatrixView<T,FortranStyle> Adjoint() const
	{ return LowerTriMatrixView<T,CStyle>::Adjoint(); }

	inline LowerTriMatrixView<RealType(T)> Real() const
	{ return LowerTriMatrixView<T,CStyle>::Real(); }

	inline LowerTriMatrixView<RealType(T)> Imag() const
	{ return LowerTriMatrixView<T,CStyle>::Imag(); }

	template <class T2> inline void operator=(const BaseMatrix<T2>&) const
	{ TMVAssert(FALSE); }

	using LowerTriMatrixView<T,CStyle>::size;

      protected :

	using GenLowerTriMatrix<T>::okij;
	using LowerTriMatrixView<T,CStyle>::ref;

    }; // FortranStyle LowerTriMatrixView

  template <class T, DiagType D, StorageType S, IndexStyle I> 
    class UpperTriMatrix : 
      public GenUpperTriMatrix<T>
    {

      public:

	//
	// Constructors
	//

#define NEW_SIZE(s) \
	itslen((s)*(s)), itsm(new T[itslen]), itss(s) \
	DEFFIRSTLAST(itsm.get(),itsm.get()+itslen)

	explicit inline UpperTriMatrix(size_t _size) : NEW_SIZE(_size) 
	{ 
	  TMVAssert(S==RowMajor || S==ColMajor); 
#ifdef TMVDEBUG
	  SetAllTo(T(888));
#endif
	}

	template <class T2> inline UpperTriMatrix(const GenMatrix<T2>& rhs) :
	  NEW_SIZE(rhs.rowsize())
	{ 
	  TMVAssert(IsReal(T2()) || IsComplex(T()));
	  TMVAssert(S==RowMajor || S==ColMajor); 
	  Copy(UpperTriMatrixViewOf(rhs,D),View()); 
	}

	template <class T2> inline UpperTriMatrix(
	    const GenUpperTriMatrix<T2>& rhs) :
	  NEW_SIZE(rhs.size())
	{ 
	  TMVAssert(IsReal(T2()) || IsComplex(T()));
	  TMVAssert(S==RowMajor || S==ColMajor); 
	  if (isunit() && !rhs.isunit()) {
	    if (rhs.size() > 0)
	      Copy(rhs.OffDiag(),OffDiag());
	  }
	  else {
	    Copy(rhs,View());
	  }
	}

	inline UpperTriMatrix(size_t _size, T x) : 
	  NEW_SIZE(_size)
	{
	  TMVAssert(S==RowMajor || S==ColMajor);
	  SetAllTo(x);
	}

	inline UpperTriMatrix(const UpperTriMatrix<T,D,S,I>& rhs) :
	  itslen(rhs.itslen), itsm(new T[itslen]), itss(rhs.itss)
	    DEFFIRSTLAST(itsm.get(),itsm.get()+itslen)
	{
	  TMVAssert(S==RowMajor || S==ColMajor);
	  memmove(itsm.get(),rhs.itsm.get(),itslen*sizeof(T));
	}

	template <DiagType D2, IndexStyle I2> inline UpperTriMatrix(
	    const UpperTriMatrix<T,D2,S,I2>& rhs) : NEW_SIZE(rhs.size())
	{
	  TMVAssert(S==RowMajor || S==ColMajor);
	  memmove(itsm.get(),rhs.cptr(),itslen*sizeof(T));
	  if (D==NonUnitDiag && D2==UnitDiag) diag().SetAllTo(T(1));
	}

	template <IndexStyle I2> inline UpperTriMatrix(
	    const Matrix<T,S,I2>& rhs) :
	  NEW_SIZE(rhs.rowsize())
	{
	  TMVAssert(S==RowMajor || S==ColMajor);
	  if (rhs.IsSquare())
	    memmove(itsm.get(),rhs.cptr(),itslen*sizeof(T));
	  else
	    Copy(UpperTriMatrixViewOf(rhs,D),View());
	}

	inline UpperTriMatrix(const GenMatrix<T>& rhs) :
	  NEW_SIZE(rhs.rowsize())
	{
	  TMVAssert(S==RowMajor || S==ColMajor);
	  Copy(UpperTriMatrixViewOf(rhs,D),View());
	}

	inline UpperTriMatrix(const GenUpperTriMatrix<RealType(T)>& rhs) :
	  NEW_SIZE(rhs.size())
	{
	  TMVAssert(S==RowMajor || S==ColMajor);
	  if (isunit() && !rhs.isunit()) {
	    if (rhs.size() > 0) OffDiag() = rhs.OffDiag();
	  }
	  else rhs.AssignToU(View());
	}

	inline UpperTriMatrix(const GenUpperTriMatrix<ComplexType(T)>& rhs) :
	  NEW_SIZE(rhs.size())
	{
	  TMVAssert(S==RowMajor || S==ColMajor);
	  TMVAssert(IsComplex(T()));
	  if (isunit() && !rhs.isunit()) {
	    if (rhs.size() > 0) OffDiag() = rhs.OffDiag();
	  }
	  else rhs.AssignToU(View());
	}

	inline UpperTriMatrix(
	    const AssignableToUpperTriMatrix<RealType(T)>& m2) :
	  NEW_SIZE(m2.size())
	{
	  TMVAssert(S==RowMajor || S==ColMajor);
	  TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
	  m2.AssignToU(View());
	}

	inline UpperTriMatrix(
	    const AssignableToUpperTriMatrix<ComplexType(T)>& m2) :
	  NEW_SIZE(m2.size())
	{
	  TMVAssert(S==RowMajor || S==ColMajor);
	  TMVAssert(IsComplex(T()));
	  TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
	  m2.AssignToU(View());
	}

#undef NEW_SIZE

	virtual inline ~UpperTriMatrix() {}

	//
	// Op=
	//

	inline UpperTriMatrix<T,D,S,I>& operator=(
	    const UpperTriMatrix<T,D,S,I>& m2)
	{ 
	  TMVAssert(size() == m2.size());
	  if (&m2 != this) memmove(itsm.get(),m2.itsm.get(),itslen*sizeof(T));
	  return *this;
	}

	template <IndexStyle I2> inline UpperTriMatrix<T,D,S,I>& operator=(
	    const UpperTriMatrix<T,D,S,I2>& m2)
	{ 
	  TMVAssert(size() == m2.size());
	  if (&m2 != this) memmove(itsm.get(),m2.itsm.get(),itslen*sizeof(T));
	  return *this;
	}

	inline UpperTriMatrix<T,D,S,I>& operator=(
	    const GenUpperTriMatrix<RealType(T)>& m2)
	{ 
	  TMVAssert(size() == m2.size());
	  TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
	  m2.AssignToU(View());
	  return *this;
	}

	inline UpperTriMatrix<T,D,S,I>& operator=(
	    const GenUpperTriMatrix<ComplexType(T)>& m2)
	{ 
	  TMVAssert(size() == m2.size());
	  TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
	  TMVAssert(IsComplex(T()));
	  m2.AssignToU(View());
	  return *this;
	}

	template <class T2> inline UpperTriMatrix<T,D,S,I>& operator=(
	    const GenUpperTriMatrix<T2>& m2)
	{ 
	  TMVAssert(size() == m2.size());
	  TMVAssert(IsReal(T2()) || sComplex(T()));
	  TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
	  Copy(m2,View());
	  return *this;
	}

	inline UpperTriMatrix<T,D,S,I>& operator=(T x) 
	{ 
	  TMVAssert(!this->isunit() || x==T(1));
	  return SetToIdentity(x); 
	}

	inline UpperTriMatrix<T,D,S,I>& operator=(
	    const AssignableToUpperTriMatrix<RealType(T)>& m2)
	{ 
	  TMVAssert(size() == m2.size());
	  TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
	  m2.AssignToU(View());
	  return *this;
	}

	inline UpperTriMatrix<T,D,S,I>& operator=(
	    const AssignableToUpperTriMatrix<ComplexType(T)>& m2)
	{ 
	  TMVAssert(size() == m2.size());
	  TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
	  TMVAssert(IsComplex(T()));
	  m2.AssignToU(View());
	  return *this;
	}

	//
	// Access
	//

	inline T operator()(size_t i, size_t j) const
	{
	  if (I == CStyle) {
	    TMVAssert(i<size());
	    TMVAssert(j<size());
	  } else {
	    TMVAssert(i>0 && i<= size());
	    TMVAssert(j>0 && j<= size());
	  }
	  if (i>j) return T(0);
	  else if (i==j && D == UnitDiag) return T(1);
	  else {
	    if (I == CStyle) {
	      TMVAssert(okij(i,j));
	      return cref(i,j);
	    } else {
	      TMVAssert(okij(i-1,j-1));
	      return cref(i-1,j-1);
	    }
	  }
	}

	inline T& operator()(size_t i, size_t j) 
	{ 
	  if (I == CStyle) {
	    TMVAssert(i<size());
	    TMVAssert(j<size());
	    TMVAssert(okij(i,j));
	    return ref(i,j);
	  } else {
	    TMVAssert(i>0 && i<= size());
	    TMVAssert(j>0 && j<= size());
	    TMVAssert(okij(i-1,j-1));
	    return ref(i-1,j-1);
	  }
	}

	inline ConstVectorView<T,I> row(size_t i, size_t j1, size_t j2) const 
	{ 
	  if (I==FortranStyle) { TMVAssert(i>0 && j1>0); }
	  const size_t ix = (I == CStyle ? i : i-1); 
	  const size_t j1x = (I == CStyle ? j1 : j1-1); 
	  TMVAssert(ix<size());
	  TMVAssert(j1x<=j2 && j2<=size());
	  TMVAssert(j1x==j2 || okij(ix,j1x));
	  return ConstVectorView<T,I>(cptr()+ix*stepi()+j1x*stepj(),
	      j2-j1x,stepj(),NonConj);
	}

	inline ConstVectorView<T,I> col(size_t j, size_t i1, size_t i2) const
	{
	  if (I==FortranStyle) { TMVAssert(j>0 && i1>0); }
	  const size_t jx = (I == CStyle ? j : j-1); 
	  const size_t i1x = (I == CStyle ? i1 : i1-1); 
	  TMVAssert(jx<size());
	  TMVAssert(i1x<=i2 && i2<=size());
	  TMVAssert(i1x==i2 || okij(i2-1,jx));
	  return ConstVectorView<T,I>(cptr()+i1x*stepi()+jx*stepj(),
	      i2-i1x,stepi(),NonConj);
	}

	inline ConstVectorView<T,I> diag() const
	{
	  TMVAssert(!isunit());
	  return ConstVectorView<T,I>(cptr(),size(),stepi()+stepj(),NonConj); 
	}

	inline ConstVectorView<T,I> diag(int i) const
	{
	  TMVAssert(i<=int(size())); 
	  TMVAssert(isunit() ? i>0 : i>=0);
	  return ConstVectorView<T,I>(cptr()+i*stepj(),
	      size()-i,stepi()+stepj(),NonConj);
	}

	inline ConstVectorView<T,I> diag(int i, size_t j1, size_t j2) const
	{
	  if (I == FortranStyle) { TMVAssert(j1 > 0); }
	  const size_t j1x = (I==CStyle ? j1 : j1-1);
	  TMVAssert(isunit() ? i>0 : i>=0);
	  TMVAssert(i<=int(size())); 
	  TMVAssert(j1x <= j2);
	  TMVAssert(j2 <= size()-i);
	  const int ds = stepi()+stepj();
	  return ConstVectorView<T,I>(cptr()+i*stepj()+j1x*ds,
	      j2-j1x,ds,NonConj);
	}

	inline VectorView<T,I> row(size_t i, size_t j1, size_t j2)
	{ 
	  if (I==FortranStyle) { TMVAssert(i>0 && j1>0); }
	  const size_t ix = (I == CStyle ? i : i-1); 
	  const size_t j1x = (I == CStyle ? j1 : j1-1); 
	  TMVAssert(ix<size());
	  TMVAssert(j1x<=j2 && j2<=size());
	  TMVAssert(j1x==j2 || okij(ix,j1x));
	  return VectorView<T,I>(ptr()+ix*stepi()+j1x*stepj(),
	      j2-j1x,stepj(),NonConj FIRSTLAST);
	}

	inline VectorView<T,I> col(size_t j, size_t i1, size_t i2)
	{
	  if (I==FortranStyle) { TMVAssert(j>0 && i1>0); }
	  const size_t jx = (I == CStyle ? j : j-1); 
	  const size_t i1x = (I == CStyle ? i1 : i1-1); 
	  TMVAssert(jx<size());
	  TMVAssert(i1x<=i2 && i2<=size());
	  TMVAssert(i1x==i2 || okij(i2-1,jx));
	  return VectorView<T,I>(ptr()+i1x*stepi()+jx*stepj(),
	    i2-i1x,stepi(),NonConj FIRSTLAST);
	}

	inline VectorView<T,I> diag()
	{
	  TMVAssert(!isunit());
	  return VectorView<T,I>(ptr(),
	      size(),stepi()+stepj(),NonConj FIRSTLAST); 
	}

	inline VectorView<T,I> diag(int i)
	{
	  TMVAssert(i<=int(size())); 
	  TMVAssert(isunit() ? i>0 : i>=0);
	  return VectorView<T,I>(ptr()+i*stepj(),
	      size()-i,stepi()+stepj(),NonConj 
	      FIRSTLAST);
	}

	inline VectorView<T,I> diag(int i, size_t j1, size_t j2)
	{
	  if (I == FortranStyle) { TMVAssert(j1>0); }
	  const size_t j1x = (I==CStyle ? j1 : j1-1);
	  TMVAssert(isunit() ? i>0 : i>=0);
	  TMVAssert(i<=int(size())); 
	  TMVAssert(j1x <= j2);
	  TMVAssert(j2 <= size()-i);
	  const int ds = stepi()+stepj();
	  return VectorView<T,I>(ptr()+i*stepj()+j1x*ds,
	      j2-j1x,ds,NonConj FIRSTLAST);
	}

	//
	// Modifying Functions
	//

	inline UpperTriMatrix<T,D,S,I>& Zero() 
	{ return SetAllTo(0); }

	inline UpperTriMatrix<T,D,S,I>& SetAllTo(T x) 
	{ VectorViewOf(itsm.get(),itslen).SetAllTo(x); return *this; }

	inline UpperTriMatrix<T,D,S,I>& Clip(RealType(T) thresh)
	{ VectorViewOf(itsm.get(),itslen).Clip(thresh); return *this; }

	inline UpperTriMatrix<T,D,S,I>& ConjugateSelf() 
	{ VectorViewOf(itsm.get(),itslen).ConjugateSelf(); return *this; }

	inline UpperTriMatrix<T,D,S,I>& InvertSelf()
	{ View().InvertSelf(); return *this; }

	inline UpperTriMatrix<T,D,S,I>& SetToIdentity(T x=T(1)) 
	{
	  TMVAssert(!isunit() || x==T(1));
	  Zero(); if (!isunit()) diag().SetAllTo(x);
	  return *this;
	}

	//
	// SubMatrix
	//

	inline ConstMatrixView<T,I> SubMatrix(
	    int i1, int i2, int j1, int j2) const
	{
	  TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,1,1));
	  const int i1x = (I==CStyle ? i1 : i1-1);
	  const int j1x = (I==CStyle ? j1 : j1-1);
	  return ConstMatrixView<T,I>(cptr()+i1x*stepi()+j1x*stepj(),
	      i2-i1x, j2-j1x,stepi(),stepj(),S,NonConj);
	}

	inline ConstMatrixView<T,I> SubMatrix(int i1, int i2, int j1, int j2,
	    int istep, int jstep) const
	{
	  StorageType newstor = S==RowMajor ?
	    jstep == 1 ? RowMajor : NoMajor :
	    istep == 1 ? ColMajor : NoMajor;
	  TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	  const int i1x = (I==CStyle ? i1 : i1-1);
	  const int j1x = (I==CStyle ? j1 : j1-1);
	  const int i2x = (I==CStyle ? i2 : i2-1+istep);
	  const int j2x = (I==CStyle ? j2 : j2-1+jstep);
	  return ConstMatrixView<T,I>(cptr()+i1x*stepi()+j1x*stepj(),
	      (i2x-i1x)/istep, (j2x-j1x)/jstep, istep*stepi(), jstep*stepj(),
	      newstor, NonConj);
	}

	inline ConstVectorView<T,I> SubVector(int i, int j,
	    int istep, int jstep, size_t size) const
	{
	  TMVAssert(View().OKSubVector(i,j,istep,jstep,size));
	  const int ix = (I==CStyle ? i : i-1);
	  const int jx = (I==CStyle ? j : j-1);
	  return ConstVectorView<T,I>(cptr()+ix*stepi()+jx*stepj(),size,
	      istep*stepi()+jstep*stepj(),NonConj);
	}

	inline ConstUpperTriMatrixView<T,I> SubTriMatrix(int i1, int i2) const
	{
	  TMVAssert(View().OKSubTriMatrix(i1,i2,1));
	  const int i1x = (I==CStyle ? i1 : i1-1);
	  return ConstUpperTriMatrixView<T,I>(cptr()+i1x*(stepi()+stepj()),
	      i2-i1x,stepi(),stepj(),D,S,NonConj);
	}

	inline ConstUpperTriMatrixView<T,I> SubTriMatrix(
	    int i1, int i2, int istep) const
	{
	  TMVAssert(View().OKSubTriMatrix(i1,i2,istep));
	  const int i1x = (I==CStyle ? i1 : i1-1);
	  const int i2x = (I==CStyle ? i2 : i2-1+istep);
	  return ConstUpperTriMatrixView<T,I>(cptr()+i1x*(stepi()+stepj()),
	      (i2x-i1x)/istep, istep*stepi(),istep*stepj(), D,
	      istep==1 ? S : NoMajor, NonConj);
	}

	inline ConstUpperTriMatrixView<T,I> OffDiag(size_t noff=1) const
	{
	  TMVAssert(size() >= noff);
	  return ConstUpperTriMatrixView<T,I>(cptr()+noff*stepj(),
	      size()-noff,stepi(),stepj(),NonUnitDiag,S,NonConj);
	}

	inline ConstUpperTriMatrixView<T,I> MakeUnitDiag() const
	{
	  return ConstUpperTriMatrixView<T,I>(cptr(),
	      size(),stepi(),stepj(),UnitDiag,S,NonConj);
	}

	inline MatrixView<T,I> SubMatrix(int i1, int i2, int j1, int j2) 
	{
	  TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,1,1));
	  const int i1x = (I==CStyle ? i1 : i1-1);
	  const int j1x = (I==CStyle ? j1 : j1-1);
	  return MatrixView<T,I>(ptr()+i1x*stepi()+j1x*stepj(),
	      i2-i1x, j2-j1x, stepi(),stepj(),S,NonConj FIRSTLAST );
	}

	inline MatrixView<T,I> SubMatrix(int i1, int i2, int j1, int j2,
	    int istep, int jstep)
	{
	  StorageType newstor = S == RowMajor ?
	    jstep == 1 ? RowMajor : NoMajor :
	    istep == 1 ? ColMajor : NoMajor;
	  TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	  const int i1x = (I==CStyle ? i1 : i1-1);
	  const int j1x = (I==CStyle ? j1 : j1-1);
	  const int i2x = (I==CStyle ? i2 : i2-1+istep);
	  const int j2x = (I==CStyle ? j2 : j2-1+jstep);
	  return MatrixView<T,I>(ptr()+i1x*stepi()+j1x*stepj(),
	      (i2x-i1x)/istep, (j2x-j1x)/jstep, istep*stepi(), jstep*stepj(),
	      newstor,NonConj FIRSTLAST );
	}

	inline VectorView<T,I> SubVector(int i, int j,
	    int istep, int jstep, size_t size)
	{
	  TMVAssert(View().OKSubVector(i,j,istep,jstep,size));
	  const int ix = (I==CStyle ? i : i-1);
	  const int jx = (I==CStyle ? j : j-1);
	  return VectorView<T,I>(ptr()+ix*stepi()+jx*stepj(),size,
	      istep*stepi()+jstep*stepj(),NonConj FIRSTLAST );
	}

	inline UpperTriMatrixView<T,I> SubTriMatrix(int i1, int i2)
	{
	  TMVAssert(View().OKSubTriMatrix(i1,i2,1));
	  const int i1x = (I==CStyle ? i1 : i1-1);
	  return UpperTriMatrixView<T,I>(ptr()+i1x*(stepi()+stepj()),
	      i2-i1x,stepi(),stepj(),D,S,NonConj FIRSTLAST);
	}

	inline UpperTriMatrixView<T,I> SubTriMatrix(int i1, int i2, int istep) 
	{
	  TMVAssert(View().OKSubTriMatrix(i1,i2,istep));
	  const int i1x = (I==CStyle ? i1 : i1-1);
	  const int i2x = (I==CStyle ? i2 : i2-1+istep);
	  return UpperTriMatrixView<T,I>(ptr()+i1x*(stepi()+stepj()),
	      (i2x-i1x)/istep,istep*stepi(),istep*stepj(),D,
	      istep==1 ? S : NoMajor,NonConj FIRSTLAST);
	}

	inline UpperTriMatrixView<T,I> OffDiag(size_t noff=1)
	{
	  TMVAssert(size() >= noff);
	  return UpperTriMatrixView<T,I>(ptr()+noff*stepj(),size()-noff,
	      stepi(),stepj(),NonUnitDiag,S,NonConj FIRSTLAST);
	}

	inline UpperTriMatrixView<T,I> MakeUnitDiag()
	{
	  return UpperTriMatrixView<T,I>(ptr(),size(),
	      stepi(),stepj(),UnitDiag,S,NonConj FIRSTLAST);
	}

	inline ConstUpperTriMatrixView<T,I> View() const
	{ 
	  return ConstUpperTriMatrixView<T,I>(cptr(),size(),
	      stepi(),stepj(),D,S,NonConj);
	}

	inline ConstLowerTriMatrixView<T,I> Transpose() const
	{ 
	  return ConstLowerTriMatrixView<T,I>(cptr(),size(),
	      stepj(),stepi(),D,TransOf(S),NonConj);
	}

	inline ConstUpperTriMatrixView<T,I> Conjugate() const
	{ 
	  return ConstUpperTriMatrixView<T,I>(cptr(),size(),
	      stepi(),stepj(),D,S,ConjOf(T,NonConj));
	}

	inline ConstLowerTriMatrixView<T,I> Adjoint() const
	{ 
	  return ConstLowerTriMatrixView<T,I>(cptr(),size(),
	      stepj(),stepi(),D,TransOf(S),ConjOf(T,NonConj));
	}

	inline UpperTriMatrixView<T,I> View() 
	{ 
	  return UpperTriMatrixView<T,I>(ptr(),size(),
	      stepi(),stepj(),D,S,NonConj FIRSTLAST);
	}

	inline LowerTriMatrixView<T,I> Transpose() 
	{ 
	  return LowerTriMatrixView<T,I>(ptr(),size(),
	      stepj(),stepi(),D,TransOf(S),NonConj FIRSTLAST);
	}

	inline UpperTriMatrixView<T,I> Conjugate() 
	{ 
	  return UpperTriMatrixView<T,I>(ptr(),size(),
	      stepi(),stepj(),D,S,ConjOf(T,NonConj) FIRSTLAST);
	}

	inline LowerTriMatrixView<T,I> Adjoint() 
	{ 
	  return LowerTriMatrixView<T,I>(ptr(),size(),
	      stepj(),stepi(),D,TransOf(S),ConjOf(T,NonConj) FIRSTLAST);
	}

	inline size_t size() const { return itss; }
	inline const T* cptr() const { return itsm.get(); }
	inline T* ptr() { return itsm.get(); }
	inline int stepi() const { return S==RowMajor ? itss : 1; }
	inline int stepj() const { return S==RowMajor ? 1 : itss; }
	inline DiagType dt() const { return D; }
	inline StorageType stor() const { return S; }
	inline ConjType ct() const { return NonConj; }
	inline bool isrm() const { return S==RowMajor; }
	inline bool iscm() const { return S==ColMajor; }
	inline bool isunit() const { return D == UnitDiag; }
	inline bool isconj() const { return false; }

      protected :

	const size_t itslen;
	auto_array<T> itsm;
	const size_t itss;

#ifdef TMVFLDEBUG
      public :
	const T*const first;
	const T*const last;
      protected :
#endif

	inline bool okij(size_t i, size_t j) const
	{
	  TMVAssert(i < size());
	  TMVAssert(j < size());
	  if (isunit()) return i<j; else return i<=j;
	}

	inline T& ref(size_t i, size_t j)
	{
	  TMVAssert(i < size());
	  TMVAssert(j < size());
	  TMVAssert(okij(i,j));
	  return *(ptr() + i*stepi() + j*stepj());
	}

	inline T cref(size_t i, size_t j) const 
	{
	  TMVAssert(i < size());
	  TMVAssert(j < size());
	  return i==j && isunit() ? T(1) : i>j ? T(0) :
	    *(cptr() + i*stepi() + j*stepj());
	}

    }; // UpperTriMatrix - RowMajor

  template <class T, DiagType D, StorageType S, IndexStyle I> 
    class LowerTriMatrix : 
      public GenLowerTriMatrix<T>
    {

      public:

	//
	// Constructors
	//

#define NEW_SIZE(s) \
	itslen((s)*(s)), itsm(new T[itslen]), itss(s) \
	DEFFIRSTLAST(itsm.get(),itsm.get()+itslen)

	explicit inline LowerTriMatrix(size_t _size) : NEW_SIZE(_size) 
	{ 
	  TMVAssert(S==RowMajor || S==ColMajor); 
#ifdef TMVDEBUG
	  SetAllTo(T(888));
#endif
	}

	template <class T2> inline LowerTriMatrix(const GenMatrix<T2>& rhs) :
	  NEW_SIZE(rhs.colsize())
	{ 
	  TMVAssert(IsReal(T2()) || IsComplex(T()));
	  TMVAssert(S==RowMajor || S==ColMajor); 
	  Copy(LowerTriMatrixViewOf(rhs,D).Transpose(),Transpose()); 
	}

	template <class T2> inline LowerTriMatrix(
	    const GenLowerTriMatrix<T2>& rhs) :
	  NEW_SIZE(rhs.size())
	{ 
	  TMVAssert(IsReal(T2()) || IsComplex(T()));
	  TMVAssert(S==RowMajor || S==ColMajor); 
	  if (isunit() && !rhs.isunit()) {
	    if (rhs.size() > 0) 
	      Copy(rhs.OffDiag().Transpose(),OffDiag().Transpose());
	  }
	  else {
	    Copy(rhs.Transpose(),Transpose());
	  }
	}

	inline LowerTriMatrix(size_t _size, T x) : NEW_SIZE(_size)
	{
	  TMVAssert(S==RowMajor || S==ColMajor);
	  SetAllTo(x);
	}

	inline LowerTriMatrix(const LowerTriMatrix<T,D,S,I>& rhs) :
	  itslen(rhs.itslen), itsm(new T[itslen]), itss(rhs.itss)
	    DEFFIRSTLAST(itsm.get(),itsm.get()+itslen)
	{
	  TMVAssert(S==RowMajor || S==ColMajor);
	  memmove(itsm.get(),rhs.itsm.get(),itslen*sizeof(T));
	}

	template <DiagType D2, IndexStyle I2> inline LowerTriMatrix(
	    const LowerTriMatrix<T,D2,S,I2>& rhs) : NEW_SIZE(rhs.size())
	{
	  TMVAssert(S==RowMajor || S==ColMajor);
	  memmove(itsm.get(),rhs.cptr(),itslen*sizeof(T));
	  if (D==NonUnitDiag && D2==UnitDiag) diag().SetAllTo(T(1));
	}

	template <IndexStyle I2> inline LowerTriMatrix(
	    const Matrix<T,S,I2>& rhs) :
	  NEW_SIZE(rhs.rowsize())
	{
	  TMVAssert(S==RowMajor || S==ColMajor);
	  if (rhs.IsSquare())
	    memmove(itsm.get(),rhs.cptr(),itslen*sizeof(T));
	  else
	    Copy(LowerTriMatrixViewOf(rhs,D).Transpose(),Transpose());
	}

	inline LowerTriMatrix(const GenMatrix<T>& rhs) :
	  NEW_SIZE(rhs.rowsize())
	{
	  TMVAssert(S==RowMajor || S==ColMajor);
	  Copy(LowerTriMatrixViewOf(rhs,D).Transpose(),Transpose());
	}

	inline LowerTriMatrix(const GenLowerTriMatrix<RealType(T)>& rhs) :
	  NEW_SIZE(rhs.size())
	{
	  TMVAssert(S==RowMajor || S==ColMajor);
	  if (isunit() && !rhs.isunit()) {
	    if (rhs.size() > 0) OffDiag() = rhs.OffDiag();
	  }
	  else rhs.AssignToL(View());
	}

	inline LowerTriMatrix(const GenLowerTriMatrix<ComplexType(T)>& rhs) :
	  NEW_SIZE(rhs.size())
	{
	  TMVAssert(S==RowMajor || S==ColMajor);
	  TMVAssert(IsComplex(T()));
	  if (isunit() && !rhs.isunit()) {
	    if (rhs.size() > 0) OffDiag() = rhs.OffDiag();
	  }
	  else rhs.AssignToL(View());
	}

	inline LowerTriMatrix(
	    const AssignableToLowerTriMatrix<RealType(T)>& m2) :
	  NEW_SIZE(m2.size())
	{
	  TMVAssert(S==RowMajor || S==ColMajor);
	  TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
	  m2.AssignToL(View());
	}

	inline LowerTriMatrix(
	    const AssignableToLowerTriMatrix<ComplexType(T)>& m2) :
	  NEW_SIZE(m2.size())
	{
	  TMVAssert(S==RowMajor || S==ColMajor);
	  TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
	  TMVAssert(IsComplex(T()));
	  m2.AssignToL(View());
	}

#undef NEW_SIZE

	virtual inline ~LowerTriMatrix() {}

	//
	// Op=
	//

	inline LowerTriMatrix<T,D,S,I>& operator=(
	    const LowerTriMatrix<T,D,S,I>& m2)
	{ 
	  TMVAssert(size() == m2.size());
	  if (&m2 != this) memmove(itsm.get(),m2.itsm.get(),itslen*sizeof(T));
	  return *this;
	}

	template <IndexStyle I2> inline LowerTriMatrix<T,D,S,I>& operator=(
	    const LowerTriMatrix<T,D,S,I2>& m2)
	{ 
	  TMVAssert(size() == m2.size());
	  if (&m2 != this) memmove(itsm.get(),m2.itsm.get(),itslen*sizeof(T));
	  return *this;
	}

	inline LowerTriMatrix<T,D,S,I>& operator=(
	    const GenLowerTriMatrix<RealType(T)>& m2)
	{ 
	  TMVAssert(size() == m2.size());
	  TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
	  m2.AssignToL(View());
	  return *this;
	}

	inline LowerTriMatrix<T,D,S,I>& operator=(
	    const GenLowerTriMatrix<ComplexType(T)>& m2)
	{ 
	  TMVAssert(size() == m2.size());
	  TMVAssert(IsComplex(T()));
	  TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
	  m2.AssignToL(View());
	  return *this;
	}

	template <class T2> inline LowerTriMatrix<T,D,S,I>& operator=(
	    const GenLowerTriMatrix<T2>& m2)
	{ 
	  TMVAssert(size() == m2.size());
	  TMVAssert(IsReal(T2()) || IsComplex(T()));
	  TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
	  Copy(m2,View());
	  return *this;
	}

	inline LowerTriMatrix<T,D,S,I>& operator=(T x) 
	{ 
	  TMVAssert(!this->isunit() || x==T(1));
	  return SetToIdentity(x); 
	}

	inline LowerTriMatrix<T,D,S,I>& operator=(
	    const AssignableToLowerTriMatrix<RealType(T)>& m2)
	{ 
	  TMVAssert(size() == m2.size());
	  TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
	  m2.AssignToL(View());
	  return *this;
	}

	inline LowerTriMatrix<T,D,S,I>& operator=(
	    const AssignableToLowerTriMatrix<ComplexType(T)>& m2)
	{ 
	  TMVAssert(size() == m2.size());
	  TMVAssert(IsComplex(T()));
	  TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
	  m2.AssignToL(View());
	  return *this;
	}

	//
	// Access
	//

	inline T operator()(size_t i, size_t j) const
	{
	  if (I == CStyle) {
	    TMVAssert(i<size());
	    TMVAssert(j<size());
	  } else {
	    TMVAssert(i>0 && i<= size());
	    TMVAssert(j>0 && j<= size());
	  }
	  if (i<j) return T(0);
	  else if (i==j && D == UnitDiag) return T(1);
	  else {
	    if (I == CStyle) {
	      TMVAssert(okij(i,j));
	      return cref(i,j);
	    } else {
	      TMVAssert(okij(i-1,j-1));
	      return cref(i-1,j-1);
	    }
	  }
	}

	inline T& operator()(size_t i, size_t j) 
	{ 
	  if (I == CStyle) {
	    TMVAssert(i<size());
	    TMVAssert(j<size());
	    TMVAssert(okij(i,j));
	    return ref(i,j);
	  } else {
	    TMVAssert(i>0 && i<= size());
	    TMVAssert(j>0 && j<= size());
	    TMVAssert(okij(i-1,j-1));
	    return ref(i-1,j-1);
	  }
	}

	inline ConstVectorView<T,I> row(size_t i, size_t j1, size_t j2) const 
	{ 
	  if (I==FortranStyle) { TMVAssert(i>0 && j1>0); }
	  const size_t ix = (I == CStyle ? i : i-1); 
	  const size_t j1x = (I == CStyle ? j1 : j1-1); 
	  TMVAssert(ix<size());
	  TMVAssert(j1x<=j2 && j2<=size());
	  TMVAssert(j1x==j2 || okij(ix,j2-1));
	  return ConstVectorView<T,I>(cptr()+ix*stepi()+j1x*stepj(),
	      j2-j1x,stepj(),NonConj);
	}

	inline ConstVectorView<T,I> col(size_t j, size_t i1, size_t i2) const
	{
	  if (I==FortranStyle) { TMVAssert(j>0 && i1>0); }
	  const size_t jx = (I == CStyle ? j : j-1); 
	  const size_t i1x = (I == CStyle ? i1 : i1-1); 
	  TMVAssert(jx<size());
	  TMVAssert(i1x<=i2 && i2<=size());
	  TMVAssert(i1x==i2 || okij(i1x,jx));
	  return ConstVectorView<T,I>(cptr()+i1x*stepi()+jx*stepj(),
	      i2-i1x,stepi(),NonConj);
	}

	inline ConstVectorView<T,I> diag() const
	{
	  TMVAssert(!isunit());
	  return ConstVectorView<T,I>(cptr(),size(),stepi()+stepj(),NonConj); 
	}

	inline ConstVectorView<T,I> diag(int i) const
	{
	  TMVAssert(-i<=int(size())); 
	  TMVAssert(isunit() ? i<0 : i<=0);
	  return ConstVectorView<T,I>(cptr()-i*stepi(),size()+i,stepi()+stepj(),
	      NonConj);
	}

	inline ConstVectorView<T,I> diag(int i, size_t j1, size_t j2) const
	{
	  if (I == FortranStyle) { TMVAssert(j1>0); }
	  const size_t j1x = (I==CStyle ? j1 : j1-1);
	  TMVAssert(isunit() ? i<0 : i<=0);
	  TMVAssert(-i<=int(size())); 
	  TMVAssert(j1x <= j2);
	  TMVAssert(j2 <= size()+i);
	  const int ds = stepi()+stepj();
	  return ConstVectorView<T,I>(cptr()-i*stepi()+j1x*ds,
	      j2-j1x,ds,NonConj);
	}

	inline VectorView<T,I> row(size_t i, size_t j1, size_t j2)
	{ 
	  if (I==FortranStyle) { TMVAssert(i>0 && j1>0); }
	  const size_t ix = (I == CStyle ? i : i-1); 
	  const size_t j1x = (I == CStyle ? j1 : j1-1); 
	  TMVAssert(ix<size());
	  TMVAssert(j1x<=j2 && j2<=size());
	  TMVAssert(j1x==j2 || okij(ix,j2-1));
	  return VectorView<T,I>(ptr()+ix*stepi()+j1x*stepj(),
	      j2-j1x,stepj(),NonConj FIRSTLAST);
	}

	inline VectorView<T,I> col(size_t j, size_t i1, size_t i2)
	{
	  if (I==FortranStyle) { TMVAssert(j>0 && i1>0); }
	  const size_t jx = (I == CStyle ? j : j-1); 
	  const size_t i1x = (I == CStyle ? i1 : i1-1); 
	  TMVAssert(jx<size());
	  TMVAssert(i1x<=i2 && i2<=size());
	  TMVAssert(i1x==i2 || okij(i1x,jx));
	  return VectorView<T,I>(ptr()+i1x*stepi()+jx*stepj(),
	      i2-i1x,stepi(),NonConj FIRSTLAST);
	}

	inline VectorView<T,I> diag()
	{
	  TMVAssert(!isunit());
	  return VectorView<T,I>(ptr(),
	      size(),stepi()+stepj(),NonConj FIRSTLAST); 
	}

	inline VectorView<T,I> diag(int i)
	{
	  TMVAssert(-i<=int(size())); 
	  TMVAssert(isunit() ? i<0 : i<=0);
	  return VectorView<T,I>(ptr()-i*stepi(),
	      size()+i,stepi()+stepj(),NonConj FIRSTLAST);
	}

	inline VectorView<T,I> diag(int i, size_t j1, size_t j2)
	{
	  if (I == FortranStyle) { TMVAssert(j1>0); }
	  const size_t j1x = (I==CStyle ? j1 : j1-1);
	  TMVAssert(isunit() ? i<0 : i<=0);
	  TMVAssert(-i<=int(size())); 
	  TMVAssert(j1x <= j2);
	  TMVAssert(j2 <= size()+i);
	  const int ds = stepi()+stepj();
	  return VectorView<T,I>(ptr()-i*stepi()+j1x*ds,
	      j2-j1x,ds,NonConj FIRSTLAST);
	}

	//
	// Modifying Functions
	//

	inline LowerTriMatrix<T,D,S,I>& Zero() 
	{ return SetAllTo(0); }

	inline LowerTriMatrix<T,D,S,I>& SetAllTo(T x) 
	{ VectorViewOf(itsm.get(),itslen).SetAllTo(x); return *this; }

	inline LowerTriMatrix<T,D,S,I>& Clip(RealType(T) thresh)
	{ VectorViewOf(itsm.get(),itslen).Clip(thresh); return *this; }

	inline LowerTriMatrix<T,D,S,I>& ConjugateSelf() 
	{ VectorViewOf(itsm.get(),itslen).ConjugateSelf(); return *this; }

	inline LowerTriMatrix<T,D,S,I>& InvertSelf()
	{ View().InvertSelf(); return *this; }

	inline LowerTriMatrix<T,D,S,I>& SetToIdentity(T x=T(1)) 
	{ 
	  TMVAssert(!isunit() || x == T(1));
	  Zero(); if (!isunit()) diag().SetAllTo(x);
	  return *this;
	}

	//
	// SubMatrix
	//

	inline ConstMatrixView<T,I> SubMatrix(
	    int i1, int i2, int j1, int j2) const
	{
	  TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,1,1));
	  const int i1x = (I==CStyle ? i1 : i1-1);
	  const int j1x = (I==CStyle ? j1 : j1-1);
	  return ConstMatrixView<T,I>(cptr()+i1x*stepi()+j1x*stepj(),
	      i2-i1x, j2-j1x,stepi(),stepj(),S,NonConj);
	}

	inline ConstMatrixView<T,I> SubMatrix(
	    int i1, int i2, int j1, int j2, int istep, int jstep) const
	{
	  StorageType newstor = S==RowMajor ?
	    jstep == 1 ? RowMajor : NoMajor :
	    istep == 1 ? ColMajor : NoMajor;
	  TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	  const int i1x = (I==CStyle ? i1 : i1-1);
	  const int j1x = (I==CStyle ? j1 : j1-1);
	  const int i2x = (I==CStyle ? i2 : i2-1+istep);
	  const int j2x = (I==CStyle ? j2 : j2-1+jstep);
	  return ConstMatrixView<T,I>(cptr()+i1x*stepi()+j1x*stepj(),
	      (i2x-i1x)/istep, (j2x-j1x)/jstep, istep*stepi(), jstep*stepj(),
	      newstor,NonConj);
	}

	inline ConstVectorView<T,I> SubVector(
	    int i, int j, int istep, int jstep, size_t size) const
	{
	  TMVAssert(View().OKSubVector(i,j,istep,jstep,size));
	  const int ix = (I==CStyle ? i : i-1);
	  const int jx = (I==CStyle ? j : j-1);
	  return ConstVectorView<T,I>(cptr()+ix*stepi()+jx*stepj(),size,
	      istep*stepi()+jstep*stepj(),NonConj);
	}

	inline ConstLowerTriMatrixView<T,I> SubTriMatrix(int i1, int i2) const
	{
	  TMVAssert(View().OKSubTriMatrix(i1,i2,1));
	  const int i1x = (I==CStyle ? i1 : i1-1);
	  return ConstLowerTriMatrixView<T,I>(cptr()+i1x*(stepi()+stepj()),
	      i2-i1x,stepi(),stepj(),D,S,NonConj);
	}

	inline ConstLowerTriMatrixView<T,I> SubTriMatrix(
	    int i1, int i2, int istep) const
	{
	  TMVAssert(View().OKSubTriMatrix(i1,i2,istep));
	  const int i1x = (I==CStyle ? i1 : i1-1);
	  const int i2x = (I==CStyle ? i2 : i2-1+istep);
	  return ConstLowerTriMatrixView<T,I>(cptr()+i1x*(stepi()+stepj()),
	      (i2x-i1x)/istep,istep*stepi(),istep*stepj(),D,
	      istep==1 ? S : NoMajor, NonConj);
	}

	inline ConstLowerTriMatrixView<T,I> OffDiag(size_t noff=1) const
	{
	  TMVAssert(size() >= noff);
	  return ConstLowerTriMatrixView<T,I>(cptr()+noff*stepi(),size()-noff,
	      stepi(),stepj(),NonUnitDiag,S,NonConj);
	}

	inline ConstLowerTriMatrixView<T,I> MakeUnitDiag() const
	{
	  return ConstLowerTriMatrixView<T,I>(cptr(),size(),
	      stepi(),stepj(),UnitDiag,S,NonConj);
	}

	inline MatrixView<T,I> SubMatrix(int i1, int i2, int j1, int j2) 
	{
	  TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,1,1));
	  const int i1x = (I==CStyle ? i1 : i1-1);
	  const int j1x = (I==CStyle ? j1 : j1-1);
	  return MatrixView<T,I>(ptr()+i1x*stepi()+j1x*stepj(),
	      i2-i1x, j2-j1x, stepi(),stepj(),S,NonConj FIRSTLAST );
	}

	inline MatrixView<T,I> SubMatrix(
	    int i1, int i2, int j1, int j2, int istep, int jstep)
	{
	  StorageType newstor = S==RowMajor ?
	    jstep == 1 ? RowMajor : NoMajor :
	    istep == 1 ? ColMajor : NoMajor;
	  TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	  const int i1x = (I==CStyle ? i1 : i1-1);
	  const int j1x = (I==CStyle ? j1 : j1-1);
	  const int i2x = (I==CStyle ? i2 : i2-1+istep);
	  const int j2x = (I==CStyle ? j2 : j2-1+jstep);
	  return MatrixView<T,I>(ptr()+i1x*stepi()+j1x*stepj(),
	      (i2x-i1x)/istep, (j2x-j1x)/jstep, istep*stepi(), jstep*stepj(),
	      newstor,NonConj FIRSTLAST );
	}

	inline VectorView<T,I> SubVector(
	    int i, int j, int istep, int jstep, size_t size)
	{
	  TMVAssert(View().OKSubVector(i,j,istep,jstep,size));
	  const int ix = (I==CStyle ? i : i-1);
	  const int jx = (I==CStyle ? j : j-1);
	  return VectorView<T,I>(ptr()+ix*stepi()+jx*stepj(),size,
	      istep*stepi()+jstep*stepj(),NonConj FIRSTLAST );
	}

	inline LowerTriMatrixView<T,I> SubTriMatrix(int i1, int i2)
	{
	  TMVAssert(View().OKSubTriMatrix(i1,i2,1));
	  const int i1x = (I==CStyle ? i1 : i1-1);
	  return LowerTriMatrixView<T,I>(ptr()+i1x*(stepi()+stepj()),
	      i2-i1x,stepi(),stepj(),D,S,NonConj FIRSTLAST);
	}

	inline LowerTriMatrixView<T,I> SubTriMatrix(int i1, int i2, int istep) 
	{
	  TMVAssert(View().OKSubTriMatrix(i1,i2,istep));
	  const int i1x = (I==CStyle ? i1 : i1-1);
	  const int i2x = (I==CStyle ? i2 : i2-1+istep);
	  return LowerTriMatrixView<T,I>(ptr()+i1x*(stepi()+stepj()),
	      (i2x-i1x)/istep,istep*stepi(),istep*stepj(),D,
	      istep==1 ? S : NoMajor,NonConj FIRSTLAST);
	}

	inline LowerTriMatrixView<T,I> OffDiag(size_t noff=1)
	{
	  TMVAssert(size() >= noff);
	  return LowerTriMatrixView<T,I>(ptr()+noff*stepi(),size()-noff,
	      stepi(),stepj(),NonUnitDiag,S,NonConj FIRSTLAST);
	}

	inline LowerTriMatrixView<T,I> MakeUnitDiag()
	{
	  return LowerTriMatrixView<T,I>(ptr(),size(),
	      stepi(),stepj(),UnitDiag,S,NonConj FIRSTLAST);
	}

	inline ConstLowerTriMatrixView<T,I> View() const
	{ 
	  return ConstLowerTriMatrixView<T,I>(cptr(),size(),
	      stepi(),stepj(),D,S,NonConj);
	}

	inline ConstUpperTriMatrixView<T,I> Transpose() const
	{ 
	  return ConstUpperTriMatrixView<T,I>(cptr(),size(),
	      stepj(),stepi(),D,TransOf(S),NonConj);
	}

	inline ConstLowerTriMatrixView<T,I> Conjugate() const
	{ 
	  return ConstLowerTriMatrixView<T,I>(cptr(),size(),
	      stepi(),stepj(),D,S,ConjOf(T,NonConj));
	}

	inline ConstUpperTriMatrixView<T,I> Adjoint() const
	{ 
	  return ConstUpperTriMatrixView<T,I>(cptr(),size(),
	      stepj(),stepi(),D,TransOf(S),ConjOf(T,NonConj));
	}

	inline LowerTriMatrixView<T,I> View() 
	{ 
	  return LowerTriMatrixView<T,I>(ptr(),size(),
	      stepi(),stepj(),D,S,NonConj FIRSTLAST);
	}

	inline UpperTriMatrixView<T,I> Transpose() 
	{ 
	  return UpperTriMatrixView<T,I>(ptr(),size(),
	      stepj(),stepi(),D,TransOf(S),NonConj FIRSTLAST);
	}

	inline LowerTriMatrixView<T,I> Conjugate() 
	{ 
	  return LowerTriMatrixView<T,I>(ptr(),size(),
	      stepi(),stepj(),D,S,ConjOf(T,NonConj) FIRSTLAST);
	}

	inline UpperTriMatrixView<T,I> Adjoint() 
	{ 
	  return UpperTriMatrixView<T,I>(ptr(),size(),
	      stepj(),stepi(),D,TransOf(S),ConjOf(T,NonConj) FIRSTLAST);
	}

	inline size_t size() const { return itss; }
	inline const T* cptr() const { return itsm.get(); }
	inline T* ptr() { return itsm.get(); }
	inline int stepi() const { return S==RowMajor ? itss : 1; }
	inline int stepj() const { return S==RowMajor ? 1 : itss; }
	inline DiagType dt() const { return D; }
	inline StorageType stor() const { return S; }
	inline ConjType ct() const { return NonConj; }
	inline bool isrm() const { return S==RowMajor; }
	inline bool iscm() const { return S==ColMajor; }
	inline bool isunit() const { return D == UnitDiag; }
	inline bool isconj() const { return false; }

      protected :

	const size_t itslen;
	auto_array<T> itsm;
	const size_t itss;

#ifdef TMVFLDEBUG
      public:
	const T*const first;
	const T*const last;
      protected :
#endif

	inline bool okij(size_t i, size_t j) const
	{
	  TMVAssert(i < size());
	  TMVAssert(j < size());
	  return isunit() ? i>j : i>=j;
	}

	inline T& ref(size_t i, size_t j)
	{
	  TMVAssert(i < size());
	  TMVAssert(j < size());
	  TMVAssert(okij(i,j));
	  return *(ptr() + i*stepi() + j*stepj());
	}

	inline T cref(size_t i, size_t j) const 
	{
	  TMVAssert(i < size());
	  TMVAssert(j < size());
	  return (i==j && isunit()) ? T(1) : (i<j) ? T(0) :
	    *(cptr() + i*stepi() + j*stepj());
	}

    }; // LowerTriMatrix - RowMajor

//---------------------------------------------------------------------------

  //
  // Special Creators: 
  //   UpperTriMatrixViewOf(m)
  //   LowerTriMatrixViewOf(m)
  //   UnitTriMatrixViewOf(t)
  //

  template <class T> inline ConstUpperTriMatrixView<T> UpperTriMatrixViewOf(
      const GenMatrix<T>& m, DiagType dt=NonUnitDiag)
  {
    TMVAssert(m.colsize()>=m.rowsize());
    return ConstUpperTriMatrixView<T>(
	m.cptr(),m.rowsize(),m.stepi(),m.stepj(),dt,m.stor(),m.ct());
  }

  template <class T, IndexStyle I> 
    inline ConstUpperTriMatrixView<T,I> UpperTriMatrixViewOf(
	const ConstMatrixView<T,I>& m, DiagType dt=NonUnitDiag)
    { 
      TMVAssert(m.colsize()>=m.rowsize());
      return ConstUpperTriMatrixView<T,I>(
	  m.cptr(),m.rowsize(),m.stepi(),m.stepj(),dt,m.stor(),m.ct());
    }

  template <class T, StorageType S, IndexStyle I>
    inline ConstUpperTriMatrixView<T,I> UpperTriMatrixViewOf(
	const Matrix<T,S,I>& m, DiagType dt=NonUnitDiag)
    {
      TMVAssert(m.colsize()>=m.rowsize());
      return ConstUpperTriMatrixView<T,I>(
	  m.cptr(),m.rowsize(),m.stepi(),m.stepj(),dt,m.stor(),m.ct());
    }

  template <class T, IndexStyle I> 
    inline UpperTriMatrixView<T,I> UpperTriMatrixViewOf(
	const MatrixView<T,I>& m, DiagType dt=NonUnitDiag)
    { 
      TMVAssert(m.colsize()>=m.rowsize());
      return UpperTriMatrixView<T,I>(m.ptr(),m.rowsize(),m.stepi(),m.stepj(),
	  dt,m.stor(),m.ct() FIRSTLAST1(m.first,m.last) );
    }

  template <class T, StorageType S, IndexStyle I>
    inline UpperTriMatrixView<T,I> UpperTriMatrixViewOf(
	Matrix<T,S,I>& m, DiagType dt=NonUnitDiag)
    {
      TMVAssert(m.colsize()>=m.rowsize());
      return UpperTriMatrixView<T,I>(m.ptr(),m.rowsize(),m.stepi(),m.stepj(),
	  dt,m.stor(),m.ct() FIRSTLAST1(m.first,m.last) );
    }

  template <class T> inline ConstUpperTriMatrixView<T> UpperTriMatrixViewOf(
      const GenUpperTriMatrix<T>& m, DiagType dt)
  {
    TMVAssert(m.colsize()>=m.rowsize());
    TMVAssert(!(m.dt()==UnitDiag && dt==NonUnitDiag));
    return ConstUpperTriMatrixView<T>(m.cptr(),m.size(),m.stepi(),m.stepj(),
	dt,m.stor(),m.ct());
  }

  template <class T, IndexStyle I> 
    inline ConstUpperTriMatrixView<T,I> UpperTriMatrixViewOf(
	const ConstUpperTriMatrixView<T,I>& m, DiagType dt)
    { 
      TMVAssert(m.colsize()>=m.rowsize());
      TMVAssert(!(m.dt()==UnitDiag && dt==NonUnitDiag));
      return ConstUpperTriMatrixView<T,I>(m.cptr(),m.size(),m.stepi(),m.stepj(),
	  dt,m.stor(),m.ct());
    }

  template <class T, DiagType D, StorageType S, IndexStyle I>
    inline ConstUpperTriMatrixView<T,I> UpperTriMatrixViewOf(
	const UpperTriMatrix<T,D,S,I>& m, DiagType dt)
    {
      TMVAssert(m.colsize()>=m.rowsize());
      TMVAssert(!(D==UnitDiag && dt==NonUnitDiag));
      return ConstUpperTriMatrixView<T,I>(m.cptr(),m.size(),m.stepi(),m.stepj(),
	  dt,m.stor(),m.ct());
    }

  template <class T, IndexStyle I> 
    inline UpperTriMatrixView<T,I> UpperTriMatrixViewOf(
	const UpperTriMatrixView<T,I>& m, DiagType dt)
    { 
      TMVAssert(m.colsize()>=m.rowsize());
      TMVAssert(!(m.dt()==UnitDiag && dt==NonUnitDiag));
      return UpperTriMatrixView<T,I>(m.ptr(),m.size(),m.stepi(),m.stepj(),
	  dt,m.stor(),m.ct() FIRSTLAST1(m.first,m.last));
    }

  template <class T, DiagType D, StorageType S, IndexStyle I>
    inline UpperTriMatrixView<T,I> UpperTriMatrixViewOf(
	UpperTriMatrix<T,D,S,I>& m, DiagType dt)
    {
      TMVAssert(m.colsize()>=m.rowsize());
      TMVAssert(!(D==UnitDiag && dt==NonUnitDiag));
      return UpperTriMatrixView<T,I>(m.ptr(),m.size(),m.stepi(),m.stepj(),
	  dt,m.stor(),m.ct() FIRSTLAST1(m.first,m.last));
    }

  template <class T> inline ConstUpperTriMatrixView<T> UpperTriMatrixViewOf(
      const T* m, size_t size, DiagType dt, StorageType stor)
  {
    TMVAssert(stor == RowMajor || stor == ColMajor);
    if (stor == RowMajor)
      return ConstUpperTriMatrixView<T>(m,size,size,1,
	  dt,RowMajor,NonConj);
    else
      return ConstUpperTriMatrixView<T>(m,size,1,size,
	  dt,ColMajor,NonConj);
  }

  template <class T> inline UpperTriMatrixView<T> UpperTriMatrixViewOf(
      T* m, size_t size, DiagType dt, StorageType stor)
  {
    TMVAssert(stor == RowMajor || stor == ColMajor);
    if (stor == RowMajor)
      return UpperTriMatrixView<T>(m,size,size,1,
	  dt,RowMajor,NonConj FIRSTLAST1(m,m+size*size));
    else
      return UpperTriMatrixView<T>(m,size,1,size,
	  dt,ColMajor,NonConj FIRSTLAST1(m,m+size*size));
  }

  template <class T> inline ConstLowerTriMatrixView<T> LowerTriMatrixViewOf(
      const GenMatrix<T>& m, DiagType dt=NonUnitDiag)
  {
    TMVAssert(m.colsize()<=m.rowsize());
    return ConstLowerTriMatrixView<T>(m.cptr(),m.colsize(),m.stepi(),m.stepj(),
	dt,m.stor(),m.ct());
  }

  template <class T, IndexStyle I> 
    inline ConstLowerTriMatrixView<T,I> LowerTriMatrixViewOf(
	const ConstMatrixView<T,I>& m, DiagType dt=NonUnitDiag)
    { 
      TMVAssert(m.colsize()<=m.rowsize());
      return ConstLowerTriMatrixView<T,I>(
	  m.cptr(),m.colsize(),m.stepi(),m.stepj(),dt,m.stor(),m.ct());
    }

  template <class T, StorageType S, IndexStyle I> 
    inline ConstLowerTriMatrixView<T,I> LowerTriMatrixViewOf(
	const Matrix<T,S,I>& m, DiagType dt=NonUnitDiag)
    {
      TMVAssert(m.colsize()<=m.rowsize());
      return ConstLowerTriMatrixView<T,I>(
	  m.cptr(),m.colsize(),m.stepi(),m.stepj(),dt,m.stor(),m.ct());
    }

  template <class T, IndexStyle I> 
    inline LowerTriMatrixView<T,I> LowerTriMatrixViewOf(
	const MatrixView<T,I>& m, DiagType dt=NonUnitDiag)
    { 
      TMVAssert(m.colsize()<=m.rowsize());
      return LowerTriMatrixView<T,I>(m.ptr(),m.colsize(),m.stepi(),m.stepj(),
	  dt,m.stor(),m.ct() FIRSTLAST1(m.first,m.last));
    }

  template <class T, StorageType S, IndexStyle I> 
    inline LowerTriMatrixView<T,I> LowerTriMatrixViewOf(
	Matrix<T,S,I>& m, DiagType dt=NonUnitDiag)
    {
      TMVAssert(m.colsize()<=m.rowsize());
      return LowerTriMatrixView<T,I>(m.ptr(),m.colsize(),m.stepi(),m.stepj(),
	  dt,m.stor(),m.ct() FIRSTLAST1(m.first,m.last));
    }

  template <class T> inline ConstLowerTriMatrixView<T> LowerTriMatrixViewOf(
      const GenLowerTriMatrix<T>& m, DiagType dt)
  {
    TMVAssert(m.colsize()<=m.rowsize());
    TMVAssert(!(m.dt()==UnitDiag && dt==NonUnitDiag));
    return ConstLowerTriMatrixView<T>(m.cptr(),m.size(),m.stepi(),m.stepj(),
	dt,m.stor(),m.ct());
  }

  template <class T, IndexStyle I> 
    inline ConstLowerTriMatrixView<T,I> LowerTriMatrixViewOf(
	const ConstLowerTriMatrixView<T,I>& m, DiagType dt)
    { 
      TMVAssert(m.colsize()<=m.rowsize());
      TMVAssert(!(m.dt()==UnitDiag && dt==NonUnitDiag));
      return ConstLowerTriMatrixView<T,I>(m.cptr(),m.size(),m.stepi(),m.stepj(),
	  dt,m.stor(),m.ct());
    }

  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline ConstLowerTriMatrixView<T,I> LowerTriMatrixViewOf(
	const LowerTriMatrix<T,D,S,I>& m, DiagType dt)
    {
      TMVAssert(m.colsize()<=m.rowsize());
      TMVAssert(!(D==UnitDiag && dt==NonUnitDiag));
      return ConstLowerTriMatrixView<T,I>(m.cptr(),m.size(),m.stepi(),m.stepj(),
	  dt,m.stor(),m.ct());
    }

  template <class T, IndexStyle I> 
    inline LowerTriMatrixView<T,I> LowerTriMatrixViewOf(
	const LowerTriMatrixView<T,I>& m, DiagType dt)
    { 
      TMVAssert(m.colsize()<=m.rowsize());
      TMVAssert(!(m.dt()==UnitDiag && dt==NonUnitDiag));
      return LowerTriMatrixView<T,I>(m.ptr(),m.size(),m.stepi(),m.stepj(),
	  dt,m.stor(),m.ct() FIRSTLAST1(m.first,m.last));
    }

  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline LowerTriMatrixView<T,I> LowerTriMatrixViewOf(
	LowerTriMatrix<T,D,S,I>& m, DiagType dt)
    {
      TMVAssert(m.colsize()<=m.rowsize());
      TMVAssert(!(D==UnitDiag && dt==NonUnitDiag));
      return LowerTriMatrixView<T,I>(m.ptr(),m.size(),m.stepi(),m.stepj(),
	  dt,m.stor(),m.ct() FIRSTLAST1(m.first,m.last));
    }

  template <class T> inline ConstLowerTriMatrixView<T> LowerTriMatrixViewOf(
      const T* m, size_t size, DiagType dt, StorageType stor)
  {
    TMVAssert(stor == RowMajor || stor == ColMajor);
    if (stor == RowMajor)
      return ConstLowerTriMatrixView<T>(m,size,size,1,
	  dt,RowMajor,NonConj);
    else
      return ConstLowerTriMatrixView<T>(m,size,1,size,
	  dt,ColMajor,NonConj);
  }

  template <class T> inline LowerTriMatrixView<T> LowerTriMatrixViewOf(
      T* m, size_t size, DiagType dt, StorageType stor)
  {
    TMVAssert(stor == RowMajor || stor == ColMajor);
    if (stor == RowMajor)
      return LowerTriMatrixView<T>(m,size,size,1,
	  dt,RowMajor,NonConj FIRSTLAST1(m,m+size*size));
    else
      return LowerTriMatrixView<T>(m,size,1,size,
	  dt,ColMajor,NonConj FIRSTLAST1(m,m+size*size));
  }

  //
  // Copy
  //

  template <class T1, class T2> inline void NonUnitDiagCopy(
      const GenUpperTriMatrix<T1>& m1, const UpperTriMatrixView<T2>& m2)
  {
    TMVAssert(IsReal(T1()) || IsComplex(T2()));
    TMVAssert(m1.size() == m2.size());
    TMVAssert(m1.dt() == NonUnitDiag);
    TMVAssert(m2.dt() == NonUnitDiag);
    const size_t N = m1.size();

    if (m1.size() > 0)
      if (m1.iscm() && m2.iscm()) 
	for(size_t j=0;j<N;++j) m2.col(j,0,j+1) = m1.col(j,0,j+1);
      else 
	for(size_t i=0;i<N;++i) m2.row(i,i,N) = m1.row(i,i,N);
  }

  template <class T1, class T2> inline void Copy(
      const GenUpperTriMatrix<T1>& m1, const UpperTriMatrixView<T2>& m2)
  {
    TMVAssert(IsReal(T1()) || IsComplex(T2()));
    TMVAssert(m1.size() == m2.size());
    TMVAssert(m1.isunit() || !m2.isunit());

    if (m1.isunit()) {
      if (m1.size() > 0)
	NonUnitDiagCopy(m1.OffDiag(),m2.OffDiag());
      if (!m2.isunit())
	m2.diag().SetAllTo(T2(1));
    } else {
      NonUnitDiagCopy(m1,m2);
    }
  }


  //
  // Swap Matrices
  //

  template <class T> void Swap(
      const UpperTriMatrixView<T>& m1, const UpperTriMatrixView<T>& m2);

  template <class T, DiagType D, StorageType S, IndexStyle I> inline void Swap(
      const UpperTriMatrixView<T>& m1, UpperTriMatrix<T,D,S,I>& m2)
  { Swap(m1,m2.View()); }

  template <class T, DiagType D, StorageType S, IndexStyle I> inline void Swap(
      UpperTriMatrix<T,D,S,I>& m1, const UpperTriMatrixView<T>& m2)
  { Swap(m1.View(),m2); }

  template <class T, DiagType D, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
    inline void Swap(
	UpperTriMatrix<T,D,S1,I1>& m1, UpperTriMatrix<T,D,S2,I2>& m2)
    { Swap(m1.View(),m2.View()); }

  template <class T> inline void Swap(
      const LowerTriMatrixView<T>& m1, const LowerTriMatrixView<T>& m2)
  { Swap(m1.Transpose(),m2.Transpose()); }

  template <class T, DiagType D, StorageType S, IndexStyle I> inline void Swap(
      const LowerTriMatrixView<T>& m1, LowerTriMatrix<T,D,S,I>& m2)
  { Swap(m1.Transpose(),m2.Transpose()); }

  template <class T, DiagType D, StorageType S, IndexStyle I> inline void Swap(
      LowerTriMatrix<T,D,S,I>& m1, const LowerTriMatrixView<T>& m2)
  { Swap(m1.Transpose(),m2.Transpose()); }

  template <class T, DiagType D, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
    inline void Swap(
	LowerTriMatrix<T,D,S1,I1>& m1, LowerTriMatrix<T,D,S2,I2>& m2)
    { Swap(m1.Transpose(),m2.Transpose()); }

  //
  // Functions of Matrices:
  //

  template <class T> inline T Det(const GenUpperTriMatrix<T>& m)
  { return m.Det(); }
  template <class T> inline T Det(const GenLowerTriMatrix<T>& m)
  { return m.Det(); }

  template <class T> inline T Trace(const GenUpperTriMatrix<T>& m)
  { return m.Trace(); }
  template <class T> inline T Trace(const GenLowerTriMatrix<T>& m)
  { return m.Trace(); }

  template <class T> inline RealType(T) Norm(const GenUpperTriMatrix<T>& m)
  { return m.Norm(); }
  template <class T> inline RealType(T) Norm(const GenLowerTriMatrix<T>& m)
  { return m.Norm(); }

  template <class T> inline RealType(T) NormF(const GenUpperTriMatrix<T>& m)
  { return m.NormF(); }
  template <class T> inline RealType(T) NormF(const GenLowerTriMatrix<T>& m)
  { return m.NormF(); }

  template <class T> inline RealType(T) Norm2(const GenUpperTriMatrix<T>& m)
  { return m.Norm2(); }
  template <class T> inline RealType(T) Norm2(const GenLowerTriMatrix<T>& m)
  { return m.Norm2(); }

  template <class T> inline RealType(T) NormInf(const GenUpperTriMatrix<T>& m)
  { return m.NormInf(); }
  template <class T> inline RealType(T) NormInf(const GenLowerTriMatrix<T>& m)
  { return m.NormInf(); }

  template <class T> inline ConstLowerTriMatrixView<T> Transpose(
      const GenUpperTriMatrix<T>& m)
  { return m.Transpose(); }
  template <class T> inline ConstUpperTriMatrixView<T> Transpose(
      const GenLowerTriMatrix<T>& m)
  { return m.Transpose(); }

  template <class T, IndexStyle I> 
    inline ConstLowerTriMatrixView<T,I> Transpose(
	const ConstUpperTriMatrixView<T,I>& m)
    { return m.Transpose(); }
  template <class T, IndexStyle I> 
    inline ConstUpperTriMatrixView<T,I> Transpose(
	const ConstLowerTriMatrixView<T,I>& m)
    { return m.Transpose(); }

  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline ConstLowerTriMatrixView<T,I> Transpose(
	const UpperTriMatrix<T,D,S,I>& m)
    { return m.Transpose(); }
  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline ConstUpperTriMatrixView<T,I> Transpose(
	const LowerTriMatrix<T,D,S,I>& m)
    { return m.Transpose(); }

  template <class T, IndexStyle I> inline LowerTriMatrixView<T,I> Transpose(
      const UpperTriMatrixView<T,I>& m)
  { return m.Transpose(); }
  template <class T, IndexStyle I> inline UpperTriMatrixView<T,I> Transpose(
      const LowerTriMatrixView<T,I>& m)
  { return m.Transpose(); }

  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline LowerTriMatrixView<T,I> Transpose(UpperTriMatrix<T,D,S,I>& m)
    { return m.Transpose(); }
  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline UpperTriMatrixView<T,I> Transpose(LowerTriMatrix<T,D,S,I>& m)
    { return m.Transpose(); }

  template <class T> inline ConstUpperTriMatrixView<T> Conjugate(
      const GenUpperTriMatrix<T>& m)
  { return m.Conjugate(); }
  template <class T> inline ConstLowerTriMatrixView<T> Conjugate(
      const GenLowerTriMatrix<T>& m)
  { return m.Conjugate(); }

  template <class T, IndexStyle I> 
    inline ConstUpperTriMatrixView<T,I> Conjugate(
	const ConstUpperTriMatrixView<T,I>& m)
    { return m.Conjugate(); }
  template <class T, IndexStyle I> 
    inline ConstLowerTriMatrixView<T,I> Conjugate(
	const ConstLowerTriMatrixView<T,I>& m)
    { return m.Conjugate(); }

  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline ConstUpperTriMatrixView<T,I> Conjugate(
	const UpperTriMatrix<T,D,S,I>& m)
    { return m.Conjugate(); }
  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline ConstLowerTriMatrixView<T,I> Conjugate(
	const LowerTriMatrix<T,D,S,I>& m)
    { return m.Conjugate(); }

  template <class T, IndexStyle I> inline UpperTriMatrixView<T,I> Conjugate(
      const UpperTriMatrixView<T,I>& m)
  { return m.Conjugate(); }
  template <class T, IndexStyle I> inline LowerTriMatrixView<T,I> Conjugate(
      const LowerTriMatrixView<T,I>& m)
  { return m.Conjugate(); }

  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline UpperTriMatrixView<T,I> Conjugate(UpperTriMatrix<T,D,S,I>& m)
    { return m.Conjugate(); }
  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline LowerTriMatrixView<T,I> Conjugate(LowerTriMatrix<T,D,S,I>& m)
    { return m.Conjugate(); }

  template <class T> inline ConstLowerTriMatrixView<T> Adjoint(
      const GenUpperTriMatrix<T>& m)
  { return m.Adjoint(); }
  template <class T> inline ConstUpperTriMatrixView<T> Adjoint(
      const GenLowerTriMatrix<T>& m)
  { return m.Adjoint(); }

  template <class T, IndexStyle I> inline LowerTriMatrixView<T,I> Adjoint(
      const ConstUpperTriMatrixView<T,I>& m)
  { return m.Adjoint(); }
  template <class T, IndexStyle I> inline UpperTriMatrixView<T,I> Adjoint(
      const ConstLowerTriMatrixView<T,I>& m)
  { return m.Adjoint(); }

  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline ConstLowerTriMatrixView<T,I> Adjoint(
	const UpperTriMatrix<T,D,S,I>& m)
    { return m.Adjoint(); }
  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline ConstUpperTriMatrixView<T,I> Adjoint(
	const LowerTriMatrix<T,D,S,I>& m)
    { return m.Adjoint(); }

  template <class T, IndexStyle I> inline LowerTriMatrixView<T,I> Adjoint(
      const UpperTriMatrixView<T,I>& m)
  { return m.Adjoint(); }
  template <class T, IndexStyle I> inline UpperTriMatrixView<T,I> Adjoint(
      const LowerTriMatrixView<T,I>& m)
  { return m.Adjoint(); }

  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline LowerTriMatrixView<T,I> Adjoint(UpperTriMatrix<T,D,S,I>& m)
    { return m.Adjoint(); }
  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline UpperTriMatrixView<T,I> Adjoint(LowerTriMatrix<T,D,S,I>& m)
    { return m.Adjoint(); }

  template <class T> inline QuotXU<T,T> Inverse(const GenUpperTriMatrix<T>& m)
  { return m.Inverse(); }
  template <class T> inline QuotXL<T,T> Inverse(const GenLowerTriMatrix<T>& m)
  { return m.Inverse(); }

  //
  // TriMatrix ==, != TriMatrix
  //

  template <class T1, class T2> bool operator==(
      const GenUpperTriMatrix<T1>& m1, const GenUpperTriMatrix<T2>& m2);
  template <class T1, class T2> inline bool operator==(
      const GenLowerTriMatrix<T1>& m1, const GenLowerTriMatrix<T2>& m2)
  { return m1.Transpose() == m2.Transpose(); }
  template <class T1, class T2> inline bool operator!=(
      const GenUpperTriMatrix<T1>& m1, const GenUpperTriMatrix<T2>& m2)
  { return !(m1 == m2); }
  template <class T1, class T2> inline bool operator!=(
      const GenLowerTriMatrix<T1>& m1, const GenLowerTriMatrix<T2>& m2)
  { return !(m1 == m2); }


  //
  // I/O
  //
 
  template <class T, DiagType D, StorageType S, IndexStyle I> 
    std::istream& operator>>(std::istream& is, 
	auto_ptr<UpperTriMatrix<T,D,S,I> >& m);
  template <class T, DiagType D, StorageType S, IndexStyle I> 
    std::istream& operator>>(std::istream& is, 
	auto_ptr<LowerTriMatrix<T,D,S,I> >& m);

  template <class T> std::istream& operator>>(
      std::istream& is, const UpperTriMatrixView<T>& m);
  template <class T> std::istream& operator>>(
      std::istream& is, const LowerTriMatrixView<T>& m);

  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline std::istream& operator>>(std::istream& is,
	UpperTriMatrix<T,D,S,I>& m)
    { return is>>m.View(); }
  template <class T, DiagType D, StorageType S, IndexStyle I> 
    inline std::istream& operator>>(std::istream& is, 
	LowerTriMatrix<T,D,S,I>& m)
    { return is>>m.View(); }

} // namespace tmv

#endif
