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
// This file defines the TMV SymMatrix and HermMatrix classes.
//
// Constructors:
//
//    SymMatrix is used for symmetric matrices, and HermMatrix is used
//    for Hermitian matrices.  For real matrices, these are the same thing:
//    A = A.Transpose().
//    But for complex, they are different:
//    A_sym = A_sym.Transpose()
//    A_herm = A_herm.Adjoint()
//
//    For these notes, I will always write SymMatrix, but (except where
//    otherwise indicated) everything applies the same for Sym and Herm.
//    Also, the Views keep track of sym/herm difference with a parameter,
//    so it is always a GenSymMatrix, ConstSymMatrixView, or 
//    SymMatrixView - never Herm in any of these.
//
//    Caveat: Complex Hermitian matrices are such that A = At, which 
//    implies that their diagonal elements are real.  Many routines
//    involving HermMatrixes assume the reality of the diagonal.
//    However, it is possible to assign a non-real value to a diagonal
//    element.  If the user does this, the results are undefined.
//
//    In addition to the type template parameter (T), SymMatrixes have two
//    additional template parameters:
//        UpLoType uplo = Upper || Lower
//        StorageType stor = RowMajor || ColMajor
//
//        They both have default values, so you can omit both, or
//        just stor.  The default values are: {Upper, RowMajor}
//
//        The first, uplo, refers to which triangular half stores the actual 
//        data.
//
//        The storage follows the same meaning as for regular Matrices.
//
//    SymMatrix<T,uplo,stor>(size_t n)
//        Makes a Symmetric Matrix with column size = row size = n
//        with _uninitialized_ values.
//
//    SymMatrix<T,uplo,stor>(size_t n, T x)
//        Makes a Symmetric Matrix with values all initialized to x
//        For Hermitian matrixces, x must be real.
//
//    SymMatrix<T,uplo,stor>(size_t n, const T* m)
//    SymMatrix<T,uplo,stor>(size_t n, const vector<T>& m)
//        Make a Symmetric Matrix with copies the elements of m which
//        fall in tha appropriate upper or lower triangle.
//        The lengths of the arrays must be n*n, but only about half the 
//        elements are used.
//
//    SymMatrix<T,uplo,stor>(const Matrix<T>& m)
//        Makes a SymMatrix which copies the corresponding elements of m.
//
//    ConstSymMatrixView<T> SymMatrixViewOf(const Matrix<T>& m, uplo)
//    ConstSymMatrixView<T> HermMatrixViewOf(const Matrix<T>& m, uplo)
//        Makes a constant SymMatrix view of the corresponding part of m.
//        While this view cannot be modified, changing the original m
//        will cause corresponding changes in this view of m.
//
//    SymMatrixView<T> SymMatrixViewOf(Matrix<T>& m, uplo)
//    SymMatrixView<T> HermMatrixViewOf(Matrix<T>& m, uplo)
//        Makes a modifiable SymMatrix view of the corresponding part of m.
//        Modifying this matrix will change the corresponding elements in
//        the original Matrix.
//
//    ConstSymMatrixView<T> SymMatrixViewOf(const T* m, size, uplo, stor)
//    ConstSymMatrixView<T> HermMatrixViewOf(const T* m, size, uplo, stor)
//    SymMatrixView<T> SymMatrixViewOf(T* m, size, uplo, stor)
//    SymMatrixView<T> HermMatrixViewOf(T* m, size, uplo, stor)
//        View the actual memory pointed to by m as a SymMatrix/HermMatrix
//        with the given size and storage.
//
// Access Functions
//
//    size_t colsize() const
//    size_t rowsize() const
//    size_t size() const
//        Return the dimensions of the SymMatrix
//
//    T& operator()(int i, int j)
//    T operator()(int i, int j) const
//        Return the (i,j) element of the SymMatrix
//
//    Vector& row(int i, int j1, int j2)
//        Return a portion of the ith row 
//        This range must be a valid range for the requested row
//        and be entirely in the upper or lower triangle.
//
//    Vector& col(int j, int i1, int i2)
//        Return a portion of the jth column
//        This range must be a valid range for the requested column
//        and be entirely in the upper or lower triangle.
//
//    Vector& diag()
//        Return the main diagonal
//
//    Vector& diag(int i, int j1, int j2)
//    Vector& diag(int i)
//        Return the super- or sub-diagonal i
//        If i > 0 return the super diagonal starting at m_0i
//        If i < 0 return the sub diagonal starting at m_|i|0
//        If j1,j2 are given, it returns the diagonal Vector 
//        either from m_j1,i+j1 to m_j2,i+j2 (for i>0) 
//        or from m_|i|+j1,j1 to m_|i|+j2,j2 (for i<0)
//
// Modifying Functions:
//
//    Zero()
//    SetAllTo(T x) 
//        For HermMatrix, x must be real.
//    Clip(RT thresh)
//    ConjugateSelf()
//    TransposeSelf()
//    SetToIdentity(x = 1)
//    SwapRowsCols(i1,i2)
//        Equivalent to swapping rows i1,i2 then swapping cols i1,i2.
//    PermuteRowsCols(const int* p)
//    ReversePermuteRowsCols(const int* p)
//        Perform a series of row/col swaps.
//    void Swap(SymMatrix& m1, SymMatrix& m2)
//        The SymMatrices must be the same size and Hermitianity.
//
// Views of a SymMatrix:
//
//    SubMatrix(int i1, int i2, int j1, int j2, int istep=1, int jstep=1)
//        This member function will return a submatrix using rows i1 to i2
//        and columns j1 to j2 which refers
//        to the same physical elements as the original.
//        The submatrix must be completely contained within the upper 
//        or lower triangle.
//
//    SubVector(int i, int j, int istep, int jstep, int size)
//        Returns a VectorView which starts at position (i,j) in the 
//        matrix, moves in the directions (istep,jstep) and has a length
//        of size.  The subvector must be completely with the upper or
//        lower triangle.
//
//    SubSymMatrix(int i1, int i2, int istep)
//        Returns the SymMatrix which runs from i1 to i2 along the diagonal
//        (not including i2) with an optional step, and includes the 
//        off diagonal in the same rows/cols.
//
//        For example, with a SymMatrix of size 10, the x's below
//        are the original data, the O's are the SubSymMatrix returned
//        with the command SubSymMatrix(3,11,2), and the #'s are the 
//        SubSymMatrix returned with SubSymMatrix(0,3)
//
//        ###xxxxxxx
//        ###xxxxxxx
//        ###xxxxxxx
//        xxxOxOxOxO
//        xxxxxxxxxx
//        xxxOxOxOxO
//        xxxxxxxxxx
//        xxxOxOxOxO
//        xxxxxxxxxx
//        xxxOxOxOxO
//
//    Transpose(m)
//    Adjoint(m)
//    Conjugate(m)
//
//
// Functions of Matrixs:
//
//    Det(m)
//    LogDet(m)
//    Trace(m)
//    Norm(m) or NormF(m)
//    NormSq(m)
//    Norm1(m) 
//    Norm2(m) 
//    NormInf(m) 
//
//    m.Inverse()
//    Inverse(m)
//    m.Inverse(minv) // Takes either a SymMatrix or Matrix argument
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
//          size 
//          ( m(0,0) )
//          ( m(1,0) m(1,1) )
//          ...
//          ( m(size-1,0) ... m(size-1,size-1) )
//
//    is >> m
//        Reads m from istream is in the compact format
//        m must already be the correct size for this to work.
//
//    is >> mptr
//        If you do not know the size of the SymMatrix to be read, you can
//        use this form where mptr is an auto_ptr to an undefined SymMatrix.
//
//
// Division Control Functions:
//
//    m.DivideUsing(dt)
//    where dt is LU, CH, or SV
//     
//    LUD(), CHD(), SVD(), and SymSVD() return the
//        corresponding Divider classes.  
//
//    For SymMatrixes, LU actually does an LDLT decomposition.
//        ie. L(Lower Triangle) * Diagonal * L.Transpose()
//        For non-positive definite matrices, D may have 2x2 blocks along
//        the diagonal, which can then be further decomposed via normal LU
//        decomposition.
//        
//    The option unique to hermitian matrixes is CH - Cholskey decomposition.
//        This is only appropriate if you know that your HermMatrix is 
//        positive definite (ie. all eigenvalues are positive).
//        (This is guaranteed, for example, if all the square 
//        submatrices have positive determinant.)
//        In this case, the SymMatrix can be decomposed into L*L.Adjoint().
//
//    Finally, a difference for the SV Decomposition for Hermitian matrices
//        is that the decomposition can be done with V = Ut.  In this case,
//        the "singular values" are really the eigenvalues of A.
//        The only caveat is that they may be negative, whereas the usual
//        definition of singular values is that they are all positive.
//        Requiring positivity would destroy V = Ut, so it seemed more
//        useful to leave them as the actual eigenvalues.  Just keep that
//        in mind if you use the singular values for anything that expects
//        them to be positive.
//
//    If m is complex, symmetric (i.e. not hermitian), then the SVD
//        does not result in an eigenvalue decomposition, and we actually
//        just do a regular SVD.  In this case, the access is through
//        SymSVD(), rather than the usual SVD() method.


#ifndef TMV_SymMatrix_H
#define TMV_SymMatrix_H

#include "tmv/TMV_BaseSymMatrix.h"
#include "tmv/TMV_BaseTriMatrix.h"
#include "tmv/TMV_BaseDiagMatrix.h"
#include "tmv/TMV_SymMatrixArithFunc.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_DiagMatrix.h"
#include <vector>

namespace tmv {

  template <class T, class T1, class T2> 
  class OProdVV;
  template <class T, class T1, class T2> 
  class ProdMM;
  template <class T, class T1, class T2> 
  class ProdUL;
  template <class T, class T1, class T2> 
  class ProdLU;

  template <class T> 
  class GenSymMatrix : 
    virtual public AssignableToSymMatrix<T>,
    public BaseMatrix<T>,
    private DivHelper<T>
  {

  public:

    //
    // Constructors
    //

    inline GenSymMatrix() {}
    inline GenSymMatrix(const GenSymMatrix<T>&) {}
    virtual inline ~GenSymMatrix() {}

    //
    // Access Functions
    //

    using AssignableToSymMatrix<T>::size;
    inline size_t colsize() const { return size(); }
    inline size_t rowsize() const { return size(); }
    using AssignableToSymMatrix<T>::sym;

    inline T operator()(int i, int j) const
    {
      TMVAssert(i>=0 && i<int(size()));
      TMVAssert(j>=0 && j<int(size()));
      return cref(i,j);
    }

    inline ConstVectorView<T> row(int i, int j1, int j2) const 
    { 
      TMVAssert(i>=0 && i<int(size()));
      TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size()));
      TMVAssert( i-j1<=0 || j2-i<=1 );
      if ((i-j1<=0 && uplo()==Upper) || (j2-i<=1 && uplo()==Lower))
        return ConstVectorView<T>(cptr()+i*stepi()+j1*stepj(),j2-j1,stepj(),
            ct()); 
      else
        return ConstVectorView<T>(cptr()+i*stepj()+j1*stepi(),j2-j1,stepi(),
            issym()?ct():ConjOf(T,ct()));
    }

    inline ConstVectorView<T> col(int j, int i1, int i2) const
    {
      TMVAssert(j>=0 && j<int(size()));
      TMVAssert(i1>=0 && i1-i2<=0 && i2<=int(size()));
      TMVAssert( j-i1<=0 || i2-j<=1 );
      if ((i2-j<=1 && uplo()==Upper) || (j-i1<=0 && uplo()==Lower))
        return ConstVectorView<T>(cptr()+i1*stepi()+j*stepj(),i2-i1,stepi(),
            ct()); 
      else 
        return ConstVectorView<T>(cptr()+i1*stepj()+j*stepi(),i2-i1,stepj(),
            issym()?ct():ConjOf(T,ct()));
    }

    inline ConstVectorView<T> diag() const
    { return ConstVectorView<T>(cptr(),size(),stepi()+stepj(),ct()); }

    inline ConstVectorView<T> diag(int i) const
    {
      TMVAssert(i>=-int(size()) && i<=int(size())); 
      if (i>=0)
        if (uplo()==Upper)
          return ConstVectorView<T>(cptr()+i*stepj(),size()-i,
              stepi()+stepj(),ct()); 
        else
          return ConstVectorView<T>(cptr()+i*stepi(),size()-i,
              stepi()+stepj(),issym()?ct():ConjOf(T,ct()));
      else
        if (uplo()==Upper)
          return ConstVectorView<T>(cptr()-i*stepj(),size()+i,
              stepi()+stepj(),issym()?ct():ConjOf(T,ct()));
        else
          return ConstVectorView<T>(cptr()-i*stepi(),size()+i,
              stepi()+stepj(),ct()); 
    }

    inline ConstVectorView<T> diag(int i, int j1, int j2) const
    {
      TMVAssert(i>=-int(size()) && i<=int(size())); 
      TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size())-std::abs(i));
      const int ds = stepi()+stepj();
      if (i>=0)
        if (uplo()==Upper)
          return ConstVectorView<T>(cptr()+i*stepj()+j1*ds,j2-j1,ds,ct()); 
        else
          return ConstVectorView<T>(cptr()+i*stepi()+j1*ds,j2-j1,ds,
              issym()?ct():ConjOf(T,ct()));
      else
        if (uplo()==Upper)
          return ConstVectorView<T>(cptr()-i*stepj()+j1*ds,j2-j1,ds,
              issym()?ct():ConjOf(T,ct()));
        else
          return ConstVectorView<T>(cptr()-i*stepi()+j1*ds,j2-j1,ds,ct()); 
    }

    template <class T2> 
    inline bool SameAs(const BaseMatrix<T2>& ) const
    { return false; }

    inline bool SameAs(const GenSymMatrix<T>& m2) const
    { 
      if (this == &m2) return true;
      else if (cptr()==m2.cptr() && size()==m2.size() && 
          (IsReal(T()) || sym() == m2.sym())) {
        if (uplo() == m2.uplo())
          return (stepi() == m2.stepi() && stepj() == m2.stepj()
              && ct() == m2.ct());
        else
          return (stepi() == m2.stepj() && stepj() == m2.stepi()
              && issym() == (ct()==m2.ct()));
      } else return false;
    }

    inline void AssignToM(const MatrixView<RealType(T)>& m2) const
    {
      TMVAssert(IsReal(T()));
      TMVAssert(m2.colsize() == size());
      TMVAssert(m2.rowsize() == size());
      AssignToS(SymMatrixViewOf(m2,Upper));
      if (size() > 0)
        m2.LowerTri().OffDiag() = 
        m2.UpperTri().OffDiag().Transpose();
    }

    inline void AssignToM(const MatrixView<ComplexType(T)>& m2) const
    {
      TMVAssert(m2.colsize() == size());
      TMVAssert(m2.rowsize() == size());
      if (issym()) {
        AssignToS(SymMatrixViewOf(m2,Upper));
        if (size() > 0)
          m2.LowerTri().OffDiag() = 
          m2.UpperTri().OffDiag().Transpose();
      } else {
        m2.diag().Imag().Zero();
        AssignToS(HermMatrixViewOf(m2,Upper));
        if (size() > 0)
          m2.LowerTri().OffDiag() = 
          m2.UpperTri().OffDiag().Adjoint();
      }
    }

    inline void AssignToS(const SymMatrixView<RealType(T)>& m2) const
    {
      TMVAssert(IsReal(T()));
      TMVAssert(m2.size() == size());
      if (!SameAs(m2)) m2.UpperTri() = UpperTri(); 
    }

    inline void AssignToS(const SymMatrixView<ComplexType(T)>& m2) const
    {
      TMVAssert(m2.size() == size());
      TMVAssert(IsReal(T()) || m2.sym() == sym());
      if (!SameAs(m2)) m2.UpperTri() = UpperTri(); 
    }

    //
    // SubMatrix
    //

    bool OKSubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) const;

    inline ConstMatrixView<T> SubMatrix(
        int i1, int i2, int j1, int j2) const
    {
      TMVAssert(OKSubMatrix(i1,i2,j1,j2,1,1));
      TMVAssert(i2-1<=j1 || j2-1<=i1);
      if ((i2-1<=j1 && uplo()==Upper) || (j2-1<=i1 && uplo()==Lower))
        return ConstMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
            i2-i1,j2-j1,stepi(),stepj(),stor(),ct());
      else
        return ConstMatrixView<T>(cptr()+i1*stepj()+j1*stepi(),
            i2-i1,j2-j1,stepj(),stepi(),TransOf(stor()),
            issym()?ct():ConjOf(T,ct()));
    }

    inline ConstMatrixView<T> SubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) const
    {
      const StorageType newstor =
      iscm() ? (istep == 1 ? ColMajor : NoMajor) :
      isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
      TMVAssert(OKSubMatrix(i1,i2,j1,j2,istep,jstep));
      TMVAssert(i2-istep<=j1 || j2-jstep<=i1);
      if ((i2-istep<=j1 && uplo()==Upper) || (j2-jstep<=i1 && uplo()==Lower))
        return ConstMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
            (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
            newstor,ct());
      else
        return ConstMatrixView<T>(cptr()+i1*stepj()+j1*stepi(),
            (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
            TransOf(newstor),issym()?ct():ConjOf(T,ct()));
    }

    bool OKSubVector(
        int i, int j, int istep, int jstep, int n) const;

    inline ConstVectorView<T> SubVector(
        int i, int j, int istep, int jstep, int n) const
    {
      TMVAssert(OKSubVector(i,j,istep,jstep,n));
      if ((i-j<=0 && uplo()==Upper) || (j-i<=0 && uplo()==Lower))
        return ConstVectorView<T>(cptr()+i*stepi()+j*stepj(),n,
            istep*stepi()+jstep*stepj(),ct());
      else
        return ConstVectorView<T>(cptr()+i*stepj()+j*stepi(),n,
            istep*stepj()+jstep*stepi(),issym()?ct():ConjOf(T,ct()));
    }

    bool OKSubSymMatrix(int i1, int i2, int istep) const;

    inline ConstSymMatrixView<T> SubSymMatrix(int i1, int i2) const
    {
      TMVAssert(OKSubSymMatrix(i1,i2,1));
      return ConstSymMatrixView<T>(cptr()+i1*(stepi()+stepj()),i2-i1,
          stepi(),stepj(),sym(),uplo(),stor(),ct());
    }

    inline ConstSymMatrixView<T> SubSymMatrix(
        int i1, int i2, int istep) const
    {
      TMVAssert(OKSubSymMatrix(i1,i2,istep));
      return ConstSymMatrixView<T>(cptr()+i1*(stepi()+stepj()),
          (i2-i1)/istep,istep*stepi(),istep*stepj(),sym(),uplo(),
          istep==1 ? stor() : NoMajor,ct());
    }

    inline ConstUpperTriMatrixView<T> UpperTri(
        DiagType dt = NonUnitDiag) const
    {
      if (uplo() == Upper)
        return ConstUpperTriMatrixView<T>(cptr(),size(),
            stepi(),stepj(),dt,stor(),ct());
      else
        return ConstUpperTriMatrixView<T>(cptr(),size(),
            stepj(),stepi(),dt,TransOf(stor()),
            issym()?ct():ConjOf(T,ct()));
    }

    inline ConstLowerTriMatrixView<T> LowerTri(
        DiagType dt = NonUnitDiag) const
    {
      if (uplo() == Lower)
        return ConstLowerTriMatrixView<T>(cptr(),size(),
            stepi(),stepj(),dt,stor(),ct());
      else
        return ConstLowerTriMatrixView<T>(cptr(),size(),
            stepj(),stepi(),dt,TransOf(stor()),
            issym()?ct():ConjOf(T,ct()));
    }

    inline ConstSymMatrixView<RealType(T)> Real() const
    {
      return ConstSymMatrixView<RealType(T)>(
          reinterpret_cast<const RealType(T)*>(cptr()),size(),
          IsReal(T()) ? stepi() : 2*stepi(),
          IsReal(T()) ? stepj() : 2*stepj(),
          Sym, uplo(), IsReal(T()) ? stor() : NoMajor,NonConj);
    }

    inline ConstSymMatrixView<RealType(T)> Imag() const
    {
      TMVAssert(IsComplex(T()));
      TMVAssert(issym());
      // The imaginary part of a Hermitian matrix is anti-symmetric
      return ConstSymMatrixView<RealType(T)>(
          reinterpret_cast<const RealType(T)*>(cptr())+1,
          size(),2*stepi(),2*stepj(), Sym,uplo(),NoMajor,NonConj);
    }

    // 
    // Views
    //

    inline ConstSymMatrixView<T> View() const
    { 
      return ConstSymMatrixView<T>(cptr(),size(),stepi(),stepj(),
          sym(),uplo(),stor(),ct());
    }

    inline ConstSymMatrixView<T> Transpose() const
    { 
      return ConstSymMatrixView<T>(cptr(),size(),stepj(),stepi(),
          sym(),UTransOf(uplo()),TransOf(stor()),ct());
    }

    inline ConstSymMatrixView<T> Conjugate() const
    { 
      return ConstSymMatrixView<T>(cptr(),size(),stepi(),stepj(),
          sym(),uplo(),stor(),ConjOf(T,ct()));
    }

    inline ConstSymMatrixView<T> Adjoint() const
    { 
      return ConstSymMatrixView<T>(cptr(),size(),stepj(),stepi(),
          sym(),UTransOf(uplo()),TransOf(stor()),ConjOf(T,ct()));
    }

    inline SymMatrixView<T> NonConst() const
    {
      return SymMatrixView<T>(const_cast<T*>(cptr()),size(),
          stepi(),stepj(),sym(),uplo(),stor(),ct()
          FIRSTLAST1(cptr(),row(colsize()-1,0,colsize()).end().GetP()));
    }

    //
    // Functions of Matrix
    //

    inline T Det() const 
    { return DivHelper<T>::Det(); }

    inline RealType(T) LogDet(T* sign=0) const
    { return DivHelper<T>::LogDet(sign); }

    inline T Trace() const
    { return diag().SumElements(); }

    inline RealType(T) Norm() const 
    { return NormF(); }

    RealType(T) NormF() const;

    RealType(T) NormSq(const RealType(T) scale = RealType(T)(1)) const;

    RealType(T) Norm1() const;

    RealType(T) DoNorm2() const;
    inline RealType(T) Norm2() const 
    {
      if (this->DivIsSet() && this->GetDivType() == SV)
        return DivHelper<T>::Norm2();
      TMV_Warning("Calling SymMatrix::Norm2 without previously calling DivideUsing(SV)");
      return DoNorm2();
    }

    inline RealType(T) NormInf() const
    { return Norm1(); }

    RealType(T) MaxAbsElement() const;

    RealType(T) DoCondition() const;
    inline RealType(T) Condition() const
    {
      if (this->DivIsSet() && this->GetDivType() == SV)
        return DivHelper<T>::Condition();
      TMV_Warning("Calling SymMatrix::Condition without previously calling DivideUsing(SV)");
      return DoCondition();
    }

    inline bool Singular() const 
    { return DivHelper<T>::Singular(); }

    template <class T1> 
    void DoInverse(const SymMatrixView<T1>& sinv) const;

    template <class T1> 
    inline void Inverse(const SymMatrixView<T1>& minv) const
    {
      TMVAssert(minv.size() == size());
      TMVAssert(isherm() == minv.isherm());
      TMVAssert(issym() == minv.issym());
      DoInverse(minv);
    }

    inline void Inverse(const MatrixView<T>& minv) const
    { DivHelper<T>::Inverse(minv); }

    template <class T1> 
    inline void Inverse(const MatrixView<T1>& minv) const
    { DivHelper<T>::Inverse(minv); }

    template <class T1, StorageType S, IndexStyle I> 
    inline void Inverse(Matrix<T1,S,I>& minv) const
    { DivHelper<T>::Inverse(minv); }

    inline void InverseATA(const MatrixView<T>& ata) const
    { DivHelper<T>::InverseATA(ata); }

    template <StorageType S, IndexStyle I> 
    inline void InverseATA(Matrix<T,S,I>& ata) const
    { DivHelper<T>::InverseATA(ata); }

    template <UpLoType U, StorageType S, IndexStyle I> 
    inline void Inverse(SymMatrix<T,U,S,I>& sinv) const
    {
      TMVAssert(issym());
      Inverse(sinv.View()); 
    }

    template <UpLoType U, StorageType S, IndexStyle I> 
    inline void Inverse(HermMatrix<T,U,S,I>& sinv) const
    {
      TMVAssert(isherm());
      Inverse(sinv.View()); 
    }

    QuotXS<T,T> QInverse() const;
    inline QuotXS<T,T> Inverse() const
    { return QInverse(); }

    auto_ptr<BaseMatrix<T> > NewCopy() const;
    auto_ptr<BaseMatrix<T> > NewView() const;
    auto_ptr<BaseMatrix<T> > NewTranspose() const;
    auto_ptr<BaseMatrix<T> > NewConjugate() const;
    auto_ptr<BaseMatrix<T> > NewAdjoint() const;
    auto_ptr<BaseMatrix<T> > NewInverse() const;

    // 
    // Division Control
    //

    using DivHelper<T>::DivideInPlace;
    using DivHelper<T>::SaveDiv;
    using DivHelper<T>::SetDiv;
    using DivHelper<T>::UnSetDiv;
    using DivHelper<T>::ReSetDiv;
    using DivHelper<T>::DivIsSet;
    using DivHelper<T>::CheckDecomp;

    inline void DivideUsing(DivType dt) const
    {
      TMVAssert(dt == CH || dt == LU || dt == SV);
      TMVAssert(isherm() || dt != CH);
      DivHelper<T>::DivideUsing(dt); 
    }

    inline const SymLDLDiv<T>& LUD() const
    {
      DivideUsing(LU);
      SetDiv();
      TMVAssert(this->GetDiv());
      TMVAssert(dynamic_cast<const SymLDLDiv<T>*>(this->GetDiv()));
      return *dynamic_cast<const SymLDLDiv<T>*>(this->GetDiv());
    }

    inline const HermCHDiv<T>& CHD() const
    {
      TMVAssert(isherm());
      DivideUsing(CH);
      SetDiv();
      TMVAssert(this->GetDiv());
      TMVAssert(dynamic_cast<const HermCHDiv<T>*>(this->GetDiv()));
      return *dynamic_cast<const HermCHDiv<T>*>(this->GetDiv());
    }

    inline const HermSVDiv<T>& SVD() const
    {
      TMVAssert(isherm());
      DivideUsing(SV);
      SetDiv();
      TMVAssert(this->GetDiv());
      TMVAssert(dynamic_cast<const HermSVDiv<T>*>(this->GetDiv()));
      return *dynamic_cast<const HermSVDiv<T>*>(this->GetDiv());
    }

    inline const SymSVDiv<T>& SymSVD() const
    {
      TMVAssert(IsComplex(T()));
      TMVAssert(issym());
      DivideUsing(SV);
      SetDiv();
      TMVAssert(this->GetDiv());
      TMVAssert(dynamic_cast<const SymSVDiv<T>*>(this->GetDiv()));
      return *dynamic_cast<const SymSVDiv<T>*>(this->GetDiv());
    }

    template <class T1> 
    inline void LDivEq(const VectorView<T1>& v) const 
    { DivHelper<T>::LDivEq(v); }
    template <class T1> 
    inline void LDivEq(const MatrixView<T1>& m) const 
    { DivHelper<T>::LDivEq(m); }
    template <class T1> 
    inline void RDivEq(const VectorView<T1>& v) const 
    { DivHelper<T>::RDivEq(v); }
    template <class T1> 
    inline void RDivEq(const MatrixView<T1>& m) const 
    { DivHelper<T>::RDivEq(m); }
    template <class T1, class T0> 
    inline void LDiv(
        const GenVector<T1>& v1, const VectorView<T0>& v0) const
    { DivHelper<T>::LDiv(v1,v0); }
    template <class T1, class T0> 
    inline void LDiv(
        const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
    { DivHelper<T>::LDiv(m1,m0); }
    template <class T1, class T0> 
    inline void RDiv(
        const GenVector<T1>& v1, const VectorView<T0>& v0) const
    { DivHelper<T>::RDiv(v1,v0); }
    template <class T1, class T0> 
    inline void RDiv(
        const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
    { DivHelper<T>::RDiv(m1,m0); }

    //
    // I/O
    //

    void Write(std::ostream& os) const;
    void WriteCompact(std::ostream& os) const;
    void Write(std::ostream& os, RealType(T) thresh) const;
    void WriteCompact(std::ostream& os, RealType(T) thresh) const;

    virtual const T* cptr() const = 0;
    virtual int stepi() const = 0;
    virtual int stepj() const = 0;
    virtual UpLoType uplo() const = 0;
    virtual StorageType stor() const = 0;
    virtual ConjType ct() const = 0;
    virtual inline bool isrm() const { return stor() == RowMajor; }
    virtual inline bool iscm() const { return stor() == ColMajor; }
    inline bool isconj() const
    { 
      TMVAssert(IsComplex(T()) || ct()==NonConj);
      return IsComplex(T()) && ct() == Conj; 
    }
    inline bool isherm() const { return IsReal(T()) || sym() == Herm; }
    inline bool issym() const { return IsReal(T()) || sym() == Sym; }
    inline bool isupper() const { return uplo() == Upper; }

    inline bool HermOK() const
    { 
      if (issym()) return true;
      else return diag().Imag().NormInf() == RealType(T)(0);
    }

    virtual T cref(int i, int j) const;

  protected :

    using DivHelper<T>::GetDiv;

    void NewDivider() const;
    inline const BaseMatrix<T>& GetMatrix() const { return *this; }

  private :

    inline GenSymMatrix<T>& operator=(const GenSymMatrix<T>&) 
    { TMVAssert(FALSE); return *this; }

  }; // GenSymMatrix

  template <class T, IndexStyle I> 
  class ConstSymMatrixView : 
    public GenSymMatrix<T>
  {
  public :

    inline ConstSymMatrixView(const ConstSymMatrixView<T,I>& rhs) :
      itsm(rhs.itsm), itss(rhs.itss), itssi(rhs.itssi), itssj(rhs.itssj),
      itssym(rhs.itssym), itsuplo(rhs.itsuplo), itsstor(rhs.itsstor),
      itsct(rhs.itsct) 
    {
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    inline ConstSymMatrixView(const GenSymMatrix<T>& rhs) :
      itsm(rhs.cptr()), itss(rhs.size()), 
      itssi(rhs.stepi()), itssj(rhs.stepj()),
      itssym(rhs.sym()), itsuplo(rhs.uplo()), itsstor(rhs.stor()),
      itsct(rhs.ct()) 
    {
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    inline ConstSymMatrixView(
        const T* _m, size_t _s, int _si, int _sj,
        SymType _sym, UpLoType _uplo, StorageType _stor,
        ConjType _ct) : 
      itsm(_m), itss(_s), itssi(_si), itssj(_sj),
      itssym(_sym), itsuplo(_uplo), itsstor(_stor), itsct(_ct)
    { 
      TMVAssert(_stor==RowMajor ? _sj == 1 : _stor==ColMajor ?
          _si==1 : true);
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    virtual inline ~ConstSymMatrixView()
    {
#ifdef TMVDEBUG
      const_cast<const T*&>(itsm) = 0;
#endif
    }

    inline size_t size() const { return itss; }
    inline const T* cptr() const { return itsm; }
    inline int stepi() const { return itssi; }
    inline int stepj() const { return itssj; }
    inline SymType sym() const { return itssym; }
    inline UpLoType uplo() const { return itsuplo; }
    inline StorageType stor() const { return itsstor; }
    inline ConjType ct() const { return itsct; }

  protected :

    const T*const itsm;
    const size_t itss;
    const int itssi;
    const int itssj;

    const SymType itssym;
    const UpLoType itsuplo;
    const StorageType itsstor;
    const ConjType itsct;

  private :

    inline ConstSymMatrixView<T>& operator=(const ConstSymMatrixView<T,I>&) 
    { TMVAssert(FALSE); return *this; }

  }; // ConstSymMatrixView

  template <class T> 
  class ConstSymMatrixView<T,FortranStyle> : 
  public ConstSymMatrixView<T,CStyle>
  {
  public :

    inline ConstSymMatrixView(const ConstSymMatrixView<T,FortranStyle>& rhs) :
      ConstSymMatrixView<T,CStyle>(rhs) {}

    inline ConstSymMatrixView(const ConstSymMatrixView<T,CStyle>& rhs) :
      ConstSymMatrixView<T,CStyle>(rhs) {}

    inline ConstSymMatrixView(const GenSymMatrix<T>& rhs) :
      ConstSymMatrixView<T,CStyle>(rhs) {}

    inline ConstSymMatrixView(const T* _m, size_t _s, int _si, int _sj,
        SymType _sym, UpLoType _uplo, StorageType _stor,
        ConjType _ct) : 
      ConstSymMatrixView<T,CStyle>(_m,_s,_si,_sj,_sym,_uplo,_stor,_ct) {}

    virtual inline ~ConstSymMatrixView() {}

    //
    // Access Functions
    //

    inline T operator()(int i, int j) const
    {
      TMVAssert(i>0 && i<=int(this->size()));
      TMVAssert(j>0 && j<=int(this->size()));
      return GenSymMatrix<T>::cref(i-1,j-1);
    }

    inline ConstVectorView<T,FortranStyle> row(
        int i, int j1, int j2) const 
    { 
      TMVAssert(i>0 && i<=int(this->size()));
      TMVAssert(j1>0 && j1-j2<=0 && j2<=int(this->size()));
      return GenSymMatrix<T>::row(i-1,j1-1,j2);
    }

    inline ConstVectorView<T,FortranStyle> col(
        int j, int i1, int i2) const
    {
      TMVAssert(j>0 && j<=int(this->size()));
      TMVAssert(i1>0 && i1-i2<=0 && i2<=int(this->size()));
      return GenSymMatrix<T>::col(j-1,i1-1,i2);
    }

    inline ConstVectorView<T,FortranStyle> diag() const
    { return GenSymMatrix<T>::diag(); }

    inline ConstVectorView<T,FortranStyle> diag(int i) const
    { return GenSymMatrix<T>::diag(i); }

    inline ConstVectorView<T,FortranStyle> diag(
        int i, int j1, int j2) const
    {
      TMVAssert(j1>0);
      return GenSymMatrix<T>::diag(i,j1-1,j2);
    }

    //
    // SubMatrix
    //

    bool OKSubMatrix(int i1, int i2, int j1, int j2, 
        int istep, int jstep) const;

    bool OKSubVector(int i, int j, int istep, int jstep, int n) const;

    bool OKSubSymMatrix(int i1, int i2, int istep) const;

    inline ConstMatrixView<T,FortranStyle> SubMatrix(
        int i1, int i2, int j1, int j2) const
    {
      TMVAssert(OKSubMatrix(i1,i2,j1,j2,1,1));
      return GenSymMatrix<T>::SubMatrix(i1-1,i2,j1-1,j2);
    }

    inline ConstMatrixView<T,FortranStyle> SubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) const
    {
      TMVAssert(OKSubMatrix(i1,i2,j1,j2,istep,jstep));
      return GenSymMatrix<T>::SubMatrix(i1-1,i2-1+istep,j1-1,j2-1+jstep,
          istep,jstep);
    }

    inline ConstVectorView<T,FortranStyle> SubVector(
        int i, int j, int istep, int jstep, int n) const
    {
      TMVAssert(OKSubVector(i,j,istep,jstep,n));
      return GenSymMatrix<T>::SubVector(i-1,j-1,istep,jstep,n);
    }

    inline ConstSymMatrixView<T,FortranStyle> SubSymMatrix(
        int i1, int i2) const
    {
      TMVAssert(OKSubSymMatrix(i1,i2,1));
      return GenSymMatrix<T>::SubSymMatrix(i1-1,i2);
    }

    inline ConstSymMatrixView<T,FortranStyle> SubSymMatrix(
        int i1, int i2, int istep) const
    {
      TMVAssert(OKSubSymMatrix(i1,i2,istep));
      return GenSymMatrix<T>::SubSymMatrix(i1-1,i2-1+istep,istep);
    }

    inline ConstUpperTriMatrixView<T,FortranStyle> UpperTri(
        DiagType dt = NonUnitDiag) const
    { return GenSymMatrix<T>::UpperTri(dt); }

    inline ConstLowerTriMatrixView<T,FortranStyle> LowerTri(
        DiagType dt = NonUnitDiag) const
    { return GenSymMatrix<T>::LowerTri(dt); }

    inline ConstSymMatrixView<RealType(T),FortranStyle> Real() const
    { return GenSymMatrix<T>::Real(); }

    inline ConstSymMatrixView<RealType(T),FortranStyle> Imag() const
    { return GenSymMatrix<T>::Imag(); }

    // 
    // Views
    //

    inline ConstSymMatrixView<T,FortranStyle> View() const
    { return GenSymMatrix<T>::View(); }

    inline ConstSymMatrixView<T,FortranStyle> Transpose() const
    { return GenSymMatrix<T>::Transpose(); }

    inline ConstSymMatrixView<T,FortranStyle> Conjugate() const
    { return GenSymMatrix<T>::Conjugate(); }

    inline ConstSymMatrixView<T,FortranStyle> Adjoint() const
    { return GenSymMatrix<T>::Adjoint(); }

    inline SymMatrixView<T,FortranStyle> NonConst() const
    { return GenSymMatrix<T>::NonConst(); }

  private :

    inline ConstSymMatrixView<T,FortranStyle>& operator=(
        const ConstSymMatrixView<T,FortranStyle>&) 
    { TMVAssert(FALSE); return *this; }

  }; // FortranStyle ConstSymMatrixView

  template <class T, IndexStyle I> 
  class SymMatrixView : 
    public GenSymMatrix<T>
  {

  public:

    //
    // Constructors
    //

    inline SymMatrixView(const SymMatrixView<T,I>& rhs) : 
      itsm(rhs.itsm), itss(rhs.itss), itssi(rhs.itssi), itssj(rhs.itssj),
      itssym(rhs.sym()), itsuplo(rhs.uplo()), itsstor(rhs.stor()),
      itsct(rhs.ct()) DEFFIRSTLAST(rhs.first,rhs.last) 
    {
#ifdef XTEST
      //TMVAssert(this->HermOK());
#endif
    }

    inline SymMatrixView(T* _m, size_t _s, int _si, int _sj, 
        SymType _sym, UpLoType _uplo, StorageType _stor, ConjType _ct 
        PARAMFIRSTLAST(T) ) :
      itsm(_m), itss(_s), itssi(_si), itssj(_sj),
      itssym(_sym), itsuplo(_uplo), itsstor(_stor), itsct(_ct)
                                                              DEFFIRSTLAST(_first,_last)
    {
      TMVAssert(_stor==RowMajor ? _sj==1 : _stor==ColMajor ?
          _si==1 : true); 
#ifdef XTEST
      //TMVAssert(this->HermOK());
#endif
    }

    virtual inline ~SymMatrixView() 
    {
#ifdef TMVDEBUG
      const_cast<T*&>(itsm) = 0;
#endif
    }

    //
    // Op=
    //

    inline const SymMatrixView<T,I>& operator=(
        const SymMatrixView<T,I>& m2) const
    { 
      TMVAssert(size() == m2.size());
      TMVAssert(IsReal(T()) || m2.sym() == sym());
      if (!SameAs(m2)) UpperTri() = m2.UpperTri(); 
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this; 
    }

    inline const SymMatrixView<T,I>& operator=(
        const SymMatrixView<T,I>& m2)
    { 
      TMVAssert(size() == m2.size());
      TMVAssert(IsReal(T()) || m2.sym() == sym());
      if (!SameAs(m2)) UpperTri() = m2.UpperTri(); 
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this; 
    }

    inline const SymMatrixView<T,I>& operator=(
        const GenSymMatrix<RealType(T)>& m2) const
    { 
      TMVAssert(size() == m2.size());
      m2.AssignToS(*this);
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this; 
    }

    inline const SymMatrixView<T,I>& operator=(
        const GenSymMatrix<ComplexType(T)>& m2) const
    { 
      TMVAssert(IsComplex(T()));
      TMVAssert(size() == m2.size());
      TMVAssert(m2.sym() == sym());
      m2.AssignToS(*this);
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this; 
    }

    template <class T2> 
    inline const SymMatrixView<T,I>& operator=(
        const GenSymMatrix<T2>& m2) const
    { 
      TMVAssert(size() == m2.size());
      TMVAssert(IsComplex(T()) || IsReal(T2()));
      TMVAssert(IsReal(T2()) || m2.sym() == sym());
      if (!SameAs(m2)) UpperTri() = m2.UpperTri(); 
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this; 
    }

    inline const SymMatrixView<T,I>& operator=(T x) const 
    { 
      TMVAssert(issym() || IMAG(x) == RealType(T)(0));
      return SetToIdentity(x); 
    }

    inline const SymMatrixView<T,I>& operator=(
        const AssignableToSymMatrix<RealType(T)>& m2) const
    { 
      TMVAssert(size() == m2.size());
      m2.AssignToS(View());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    inline const SymMatrixView<T,I>& operator=(
        const AssignableToSymMatrix<ComplexType(T)>& m2) const
    { 
      TMVAssert(IsComplex(T()));
      TMVAssert(size() == m2.size());
      TMVAssert(sym() == m2.sym());
      m2.AssignToS(View());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    inline const SymMatrixView<T,I>& operator=(
        const GenDiagMatrix<RealType(T)>& m2) const
    { 
      TMVAssert(size() == m2.size());
      m2.AssignToD(DiagMatrixViewOf(diag()));
      UpperTri().OffDiag().Zero();
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    inline const SymMatrixView<T,I>& operator=(
        const GenDiagMatrix<ComplexType(T)>& m2) const
    { 
      TMVAssert(IsComplex(T()));
      TMVAssert(size() == m2.size());
      TMVAssert(issym());
      m2.AssignToD(DiagMatrixViewOf(diag()));
      UpperTri().OffDiag().Zero();
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    template <class T2> 
    inline const SymMatrixView<T,I>& operator=(
        const OProdVV<RealType(T),T2,T2>& opvv) const
    {
      TMVAssert(size() == opvv.colsize());
      TMVAssert(size() == opvv.rowsize());
      TMVAssert(opvv.GetV1().SameAs(issym() ? opvv.GetV2().View() :
            opvv.GetV2().Conjugate()));
      TMVAssert(issym() || IMAG(opvv.GetX()) == RealType(T)(0));
      Rank1Update<false>(opvv.GetX(),opvv.GetV1(),*this);
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    template <class T2> 
    inline const SymMatrixView<T,I>& operator=(
        const OProdVV<ComplexType(T),T2,T2>& opvv) const
    {
      TMVAssert(IsComplex(T()));
      TMVAssert(size() == opvv.colsize());
      TMVAssert(size() == opvv.rowsize());
      TMVAssert(opvv.GetV1().SameAs(issym() ? opvv.GetV2().View() :
            opvv.GetV2().Conjugate()));
      TMVAssert(issym() || IMAG(opvv.GetX()) == RealType(T)(0));
      Rank1Update<false>(opvv.GetX(),opvv.GetV1(),*this);
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    template <class T2> 
    inline const SymMatrixView<T,I>& operator=(
        const ProdMM<RealType(T),T2,T2>& pmm) const
    {
      TMVAssert(size() == pmm.colsize());
      TMVAssert(size() == pmm.rowsize());
      TMVAssert(pmm.GetM1().SameAs(issym() ? pmm.GetM2().Transpose() :
            pmm.GetM2().Adjoint()));
      TMVAssert(issym() || IMAG(pmm.GetX()) == RealType(T)(0));
      RankKUpdate<false>(pmm.GetX(),pmm.GetM1(),*this);
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    template <class T2> 
    inline const SymMatrixView<T,I>& operator=(
        const ProdMM<ComplexType(T),T2,T2>& pmm) const
    {
      TMVAssert(IsComplex(T()));
      TMVAssert(size() == pmm.colsize());
      TMVAssert(size() == pmm.rowsize());
      TMVAssert(pmm.GetM1().SameAs(issym() ? pmm.GetM2().Transpose() :
            pmm.GetM2().Adjoint()));
      TMVAssert(issym() || IMAG(pmm.GetX()) == RealType(T)(0));
      RankKUpdate<false>(pmm.GetX(),pmm.GetM1(),*this);
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    template <class T2> 
    inline const SymMatrixView<T,I>& operator=(
        const ProdUL<RealType(T),T2,T2>& pmm) const
    {
      TMVAssert(size() == pmm.colsize());
      TMVAssert(size() == pmm.rowsize());
      TMVAssert(pmm.GetM1().SameAs(issym() ? pmm.GetM2().Transpose() :
            pmm.GetM2().Adjoint()));
      TMVAssert(issym() || IMAG(pmm.GetX()) == RealType(T)(0));
      RankKUpdate<false>(pmm.GetX(),pmm.GetM1(),*this);
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    template <class T2> 
    inline const SymMatrixView<T,I>& operator=(
        const ProdUL<ComplexType(T),T2,T2>& pmm) const
    {
      TMVAssert(IsComplex(T()));
      TMVAssert(size() == pmm.colsize());
      TMVAssert(size() == pmm.rowsize());
      TMVAssert(pmm.GetM1().SameAs(issym() ? pmm.GetM2().Transpose() :
            pmm.GetM2().Adjoint()));
      TMVAssert(issym() || IMAG(pmm.GetX()) == RealType(T)(0));
      RankKUpdate<false>(pmm.GetX(),pmm.GetM1(),*this);
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    template <class T2> 
    inline const SymMatrixView<T,I>& operator=(
        const ProdLU<RealType(T),T2,T2>& pmm) const
    {
      TMVAssert(size() == pmm.colsize());
      TMVAssert(size() == pmm.rowsize());
      TMVAssert(pmm.GetM1().SameAs(issym() ? pmm.GetM2().Transpose() :
            pmm.GetM2().Adjoint()));
      TMVAssert(issym() || IMAG(pmm.GetX()) == RealType(T)(0));
      RankKUpdate<false>(pmm.GetX(),pmm.GetM1(),*this);
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    template <class T2> 
    inline const SymMatrixView<T,I>& operator=(
        const ProdLU<ComplexType(T),T2,T2>& pmm) const
    {
      TMVAssert(IsComplex(T()));
      TMVAssert(size() == pmm.colsize());
      TMVAssert(size() == pmm.rowsize());
      TMVAssert(pmm.GetM1().SameAs(issym() ? pmm.GetM2().Transpose() :
            pmm.GetM2().Adjoint()));
      TMVAssert(issym() || IMAG(pmm.GetX()) == RealType(T)(0));
      RankKUpdate<false>(pmm.GetX(),pmm.GetM1(),*this);
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }


    //
    // Access
    //

    typedef RefType(T) reference;

    inline reference operator()(int i,int j) const 
    {
      TMVAssert(i>=0 && i<int(size()));
      TMVAssert(j>=0 && j<int(size()));
      return ref(i,j); 
    }

    inline VectorView<T> row(int i, int j1, int j2) const 
    { 
      TMVAssert(i>=0 && i<int(size()));
      TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size()));
      TMVAssert( i-j1<=0 || j2-i<=1 );
      if ((i-j1<=0 && uplo()==Upper) || (j2-i<=1 && uplo()==Lower))
        return VectorView<T>(ptr()+i*stepi()+j1*stepj(),j2-j1,stepj(),
            ct() FIRSTLAST); 
      else
        return VectorView<T>(ptr()+i*stepj()+j1*stepi(),j2-j1,stepi(),
            this->issym()?ct():ConjOf(T,ct()) FIRSTLAST);
    }

    inline VectorView<T> col(int j, int i1, int i2) const
    {
      TMVAssert(j>=0 && j<int(size()));
      TMVAssert(i1>=0 && i1-i2<=0 && i2<=int(size()));
      TMVAssert( j-i1<=0 || i2-j<=1 );
      if ((i2-j<=1 && uplo()==Upper) || (j-i1<=0 && uplo()==Lower))
        return VectorView<T>(ptr()+i1*stepi()+j*stepj(),i2-i1,stepi(),
            ct() FIRSTLAST); 
      else 
        return VectorView<T>(ptr()+i1*stepj()+j*stepi(),i2-i1,stepj(),
            this->issym()?ct():ConjOf(T,ct()) FIRSTLAST);
    }

    inline VectorView<T> diag() const
    { return VectorView<T>(ptr(),size(),stepi()+stepj(),ct() FIRSTLAST); }

    inline VectorView<T> diag(int i) const
    {
      TMVAssert(i>=-int(size()) && i<=int(size())); 
      if (i>=0)
        if (uplo()==Upper)
          return VectorView<T>(ptr()+i*stepj(),size()-i,
              stepi()+stepj(),ct() FIRSTLAST); 
        else
          return VectorView<T>(ptr()+i*stepi(),size()-i,
              stepi()+stepj(),this->issym()?ct():ConjOf(T,ct()) FIRSTLAST);
      else
        if (uplo()==Upper)
          return VectorView<T>(ptr()-i*stepj(),size()+i,
              stepi()+stepj(),this->issym()?ct():ConjOf(T,ct()) FIRSTLAST);
        else
          return VectorView<T>(ptr()-i*stepi(),size()+i,
              stepi()+stepj(),ct() FIRSTLAST); 
    }

    inline VectorView<T> diag(int i, int j1, int j2) const
    {
      TMVAssert(i>=-int(size()) && i<=int(size())); 
      TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size())-std::abs(i));
      const int ds = stepi()+stepj();
      if (i>=0)
        if (uplo()==Upper)
          return VectorView<T>(ptr()+i*stepj()+j1*ds,j2-j1,ds,ct() 
              FIRSTLAST); 
        else
          return VectorView<T>(ptr()+i*stepi()+j1*ds,j2-j1,ds, 
              this->issym()?ct():ConjOf(T,ct()) FIRSTLAST);
      else
        if (uplo()==Upper)
          return VectorView<T>(ptr()-i*stepj()+j1*ds,j2-j1,ds, 
              this->issym()?ct():ConjOf(T,ct()) FIRSTLAST);
        else
          return VectorView<T>(ptr()-i*stepi()+j1*ds,j2-j1,ds,ct() 
              FIRSTLAST); 
    }

    //
    // Modifying Functions
    //

    inline const SymMatrixView<T,I>& Zero() const 
    { UpperTri().Zero(); return *this; }

    inline const SymMatrixView<T,I>& SetAllTo(T x) const
    { 
      TMVAssert(IMAG(x)==RealType(T)(0) || this->issym());
      UpperTri().SetAllTo(x); return *this; 
    }

    inline const SymMatrixView<T,I>& Clip(RealType(T) thresh) const
    { UpperTri().Clip(thresh); return *this; }

    inline const SymMatrixView<T,I>& ConjugateSelf() const
    { if (IsComplex(T())) UpperTri().ConjugateSelf(); return *this; }

    inline const SymMatrixView<T,I>& TransposeSelf() const
    { if (!this->issym()) UpperTri().ConjugateSelf(); return *this; }

    inline const SymMatrixView<T,I>& SetToIdentity(T x=T(1)) const
    { 
      TMVAssert(IMAG(x)==RealType(T)(0) || this->issym());
      Zero(); diag().SetAllTo(x); return *this; 
    }

    const SymMatrixView<T,I>& SwapRowsCols(int i1, int i2) const;

    const SymMatrixView<T,I>& PermuteRowsCols(const int* p,
        int i1, int i2) const;

    const SymMatrixView<T,I>& ReversePermuteRowsCols(const int* p,
        int i1, int i2) const;

    inline const SymMatrixView<T,I>& PermuteRowsCols(const int* p) const
    { return PermuteRowsCols(p,0,size()); }

    inline const SymMatrixView<T,I>& ReversePermuteRowsCols(
        const int* p) const
    { return ReversePermuteRowsCols(p,0,size()); }


    //
    // SubMatrix
    //

    inline MatrixView<T> SubMatrix(int i1, int i2, int j1, int j2) const
    {
      TMVAssert(GenSymMatrix<T>::OKSubMatrix(i1,i2,j1,j2,1,1));
      TMVAssert(i2-1<=j1 || j2-1<=i1);
      if ((i2-j1<=1 && uplo()==Upper) || (j2-i1<=1 && uplo()==Lower))
        return MatrixView<T>(ptr()+i1*stepi()+j1*stepj(),
            i2-i1,j2-j1,stepi(),stepj(),stor(),ct() FIRSTLAST);
      else
        return MatrixView<T>(ptr()+i1*stepj()+j1*stepi(),
            i2-i1,j2-j1,stepj(),stepi(),TransOf(stor()),
            this->issym()?ct():ConjOf(T,ct()) FIRSTLAST);
    }

    inline MatrixView<T> SubMatrix(int i1, int i2, int j1, int j2,
        int istep, int jstep) const
    {
      const StorageType newstor =
      this->iscm() ? (istep == 1 ? ColMajor : NoMajor) :
      this->isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
      TMVAssert(GenSymMatrix<T>::OKSubMatrix(i1,i2,j1,j2,istep,jstep));
      TMVAssert(i2-istep<=j1 || j2-jstep<=i1);
      if ((i2-istep<=j1 && uplo()==Upper) || (j2-jstep<=i1 && uplo()==Lower))
        return MatrixView<T>(ptr()+i1*stepi()+j1*stepj(),
            (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
            newstor,ct() FIRSTLAST);
      else
        return MatrixView<T>(ptr()+i1*stepj()+j1*stepi(),
            (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
            TransOf(newstor),this->issym()?ct():ConjOf(T,ct()) FIRSTLAST);
    }

    inline VectorView<T> SubVector(int i, int j,
        int istep, int jstep, int n) const
    {
      TMVAssert(GenSymMatrix<T>::OKSubVector(i,j,istep,jstep,n));
      if ((i-j<=0 && uplo()==Upper) || (j-i<=0 && uplo()==Lower))
        return VectorView<T>(ptr()+i*stepi()+j*stepj(),n,
            istep*stepi()+jstep*stepj(),ct() FIRSTLAST);
      else
        return VectorView<T>(ptr()+i*stepj()+j*stepi(),n,
            istep*stepj()+jstep*stepi(),this->issym()?ct():ConjOf(T,ct()) 
            FIRSTLAST);
    }

    inline SymMatrixView<T,I> SubSymMatrix(int i1, int i2) const
    {
      TMVAssert(GenSymMatrix<T>::OKSubSymMatrix(i1,i2,1));
      return SymMatrixView<T,I>(ptr()+i1*(stepi()+stepj()),i2-i1,
          stepi(),stepj(),sym(),uplo(),stor(),ct() FIRSTLAST);
    }

    inline SymMatrixView<T,I> SubSymMatrix(int i1, int i2,
        int istep) const
    {
      TMVAssert(GenSymMatrix<T>::OKSubSymMatrix(i1,i2,istep));
      return SymMatrixView<T,I>(ptr()+i1*(stepi()+stepj()),(i2-i1)/istep,
          istep*stepi(),istep*stepj(),sym(),uplo(),
          istep==1 ? stor() : NoMajor,ct() FIRSTLAST);
    }

    inline UpperTriMatrixView<T> UpperTri(DiagType dt = NonUnitDiag) const
    {
      if (uplo() == Upper)
        return UpperTriMatrixView<T>(ptr(),size(),
            stepi(),stepj(),dt,stor(),ct() FIRSTLAST);
      else
        return UpperTriMatrixView<T>(ptr(),size(),
            stepj(),stepi(),dt,TransOf(stor()),
            this->issym()?ct():ConjOf(T,ct()) FIRSTLAST);
    }

    inline LowerTriMatrixView<T> LowerTri(DiagType dt = NonUnitDiag) const
    {
      if (uplo() == Lower)
        return LowerTriMatrixView<T>(ptr(),size(),
            stepi(),stepj(),dt,stor(),ct() FIRSTLAST);
      else
        return LowerTriMatrixView<T>(ptr(),size(),
            stepj(),stepi(),dt,TransOf(stor()),
            this->issym()?ct():ConjOf(T,ct()) FIRSTLAST);
    }

    inline SymMatrixView<T,I> View() const
    { return *this; }

    inline SymMatrixView<T,I> Transpose() const
    {
      return SymMatrixView<T,I>(ptr(),size(),stepj(),stepi(),
          sym(),UTransOf(uplo()),TransOf(stor()),ct() FIRSTLAST);
    }

    inline SymMatrixView<T,I> Conjugate() const
    {
      return SymMatrixView<T,I>(ptr(),size(),stepi(),stepj(),
          sym(),uplo(),stor(),ConjOf(T,ct()) FIRSTLAST);
    }

    inline SymMatrixView<T,I> Adjoint() const
    {
      return SymMatrixView<T,I>(ptr(),size(),stepj(),stepi(),
          sym(),UTransOf(uplo()),TransOf(stor()),ConjOf(T,ct()) 
          FIRSTLAST);
    }

    inline SymMatrixView<RealType(T),I> Real() const
    {
      return SymMatrixView<RealType(T)>(
          reinterpret_cast<RealType(T)*>(ptr()),size(),
          IsReal(T()) ? stepi() : 2*stepi(),
          IsReal(T()) ? stepj() : 2*stepj(),
          sym(),uplo(), IsReal(T()) ? stor() : NoMajor, NonConj
#ifdef TMVFLDEBUG
          ,reinterpret_cast<const RealType(T)*>(first)
          ,reinterpret_cast<const RealType(T)*>(last)
#endif
          );
    }

    inline SymMatrixView<RealType(T),I> Imag() const
    {
      TMVAssert(IsComplex(T()));
      TMVAssert(this->issym());
      return SymMatrixView<RealType(T)>(
          reinterpret_cast<RealType(T)*>(ptr())+1,
          size(),2*stepi(),2*stepj(),sym(),uplo(),NoMajor, NonConj
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
    inline SymType sym() const { return itssym; }
    inline UpLoType uplo() const { return itsuplo; }
    inline StorageType stor() const { return itsstor; }
    inline ConjType ct() const { return itsct; }
    using GenSymMatrix<T>::issym;
    using GenSymMatrix<T>::isherm;

    reference ref(int i, int j) const;

  protected :

    T*const itsm;
    const size_t itss;
    const int itssi;
    const int itssj;

    const SymType itssym;
    const UpLoType itsuplo;
    const StorageType itsstor;
    const ConjType itsct;

#ifdef TMVFLDEBUG
  public:
    const T*const first;
    const T*const last;
  protected:
#endif

  }; // SymMatrixView

  template <class T> 
  class SymMatrixView<T,FortranStyle> : 
  public SymMatrixView<T,CStyle>
  {

  public:

    //
    // Constructors
    //

    inline SymMatrixView(const SymMatrixView<T,FortranStyle>& rhs) : 
      SymMatrixView<T,CStyle>(rhs) {}

    inline SymMatrixView(const SymMatrixView<T,CStyle>& rhs) : 
      SymMatrixView<T,CStyle>(rhs) {}

    inline SymMatrixView(
        T* _m, size_t _s, int _si, int _sj, 
        SymType _sym, UpLoType _uplo, StorageType _stor, ConjType _ct 
        PARAMFIRSTLAST(T) ) :
      SymMatrixView<T,CStyle>(_m,_s,_si,_sj,_sym,_uplo,_stor,_ct
          FIRSTLAST1(_first,_last) ) {}

    virtual inline ~SymMatrixView() {} 

    //
    // Op=
    //

    inline const SymMatrixView<T,FortranStyle>& operator=(
        const SymMatrixView<T,FortranStyle>& m2) const
    { SymMatrixView<T,CStyle>::operator=(m2); return *this; }

    inline const SymMatrixView<T,FortranStyle>& operator=(
        const SymMatrixView<T,FortranStyle>& m2)
    { SymMatrixView<T,CStyle>::operator=(m2); return *this; }

    inline const SymMatrixView<T,FortranStyle>& operator=(
        const SymMatrixView<T,CStyle>& m2) const
    { SymMatrixView<T,CStyle>::operator=(m2); return *this; }

    inline const SymMatrixView<T,FortranStyle>& operator=(
        const SymMatrixView<T,CStyle>& m2)
    { SymMatrixView<T,CStyle>::operator=(m2); return *this; }

    inline const SymMatrixView<T,FortranStyle>& operator=(
        const GenSymMatrix<RealType(T)>& m2) const
    { SymMatrixView<T,CStyle>::operator=(m2); return *this; }

    inline const SymMatrixView<T,FortranStyle>& operator=(
        const GenSymMatrix<ComplexType(T)>& m2) const
    { SymMatrixView<T,CStyle>::operator=(m2); return *this; }

    template <class T2> 
    inline const SymMatrixView<T,FortranStyle>& operator=(
        const GenSymMatrix<T2>& m2) const
    { SymMatrixView<T,CStyle>::operator=(m2); return *this; }

    inline const SymMatrixView<T,FortranStyle>& operator=(T x) const 
    { SymMatrixView<T,CStyle>::operator=(x); return *this; }

    inline const SymMatrixView<T,FortranStyle>& operator=(
        const AssignableToSymMatrix<RealType(T)>& m2) const
    { SymMatrixView<T,CStyle>::operator=(m2); return *this; }

    inline const SymMatrixView<T,FortranStyle>& operator=(
        const AssignableToSymMatrix<ComplexType(T)>& m2) const
    { SymMatrixView<T,CStyle>::operator=(m2); return *this; }

    inline const SymMatrixView<T,FortranStyle>& operator=(
        const GenDiagMatrix<RealType(T)>& m2) const
    { SymMatrixView<T,CStyle>::operator=(m2); return *this; }

    inline const SymMatrixView<T,FortranStyle>& operator=(
        const GenDiagMatrix<ComplexType(T)>& m2) const
    { SymMatrixView<T,CStyle>::operator=(m2); return *this; }

    template <class T2>
    inline const SymMatrixView<T,FortranStyle>& operator=(
        const OProdVV<T,T2,T2>& opvv) const
    { SymMatrixView<T,CStyle>::operator=(opvv); return *this; }

    template <class T2>
    inline const SymMatrixView<T,FortranStyle>& operator=(
        const ProdMM<T,T2,T2>& pmm) const
    { SymMatrixView<T,CStyle>::operator=(pmm); return *this; }

    template <class T2>
    inline const SymMatrixView<T,FortranStyle>& operator=(
        const ProdLU<T,T2,T2>& pmm) const
    { SymMatrixView<T,CStyle>::operator=(pmm); return *this; }

    template <class T2>
    inline const SymMatrixView<T,FortranStyle>& operator=(
        const ProdUL<T,T2,T2>& pmm) const
    { SymMatrixView<T,CStyle>::operator=(pmm); return *this; }


    //
    // Access
    //

    inline RefType(T) operator()(int i,int j) const 
    {
      TMVAssert(i>0 && i<=int(this->size()));
      TMVAssert(j>0 && j<=int(this->size()));
      return SymMatrixView<T,CStyle>::ref(i-1,j-1); 
    }

    inline VectorView<T,FortranStyle> row(int i, int j1, int j2) const 
    { 
      TMVAssert(i>0 && i<=int(this->size()));
      TMVAssert(j1>0 && j1-j2<=0 && j2<=int(this->size()));
      return SymMatrixView<T,CStyle>::row(i-1,j1-1,j2);
    }

    inline VectorView<T,FortranStyle> col(int j, int i1, int i2) const
    {
      TMVAssert(j>0 && j<=int(this->size()));
      TMVAssert(i1>0 && i1-i2<=0 && i2<=int(this->size()));
      return SymMatrixView<T,CStyle>::col(j-1,i1-1,i2);
    }

    inline VectorView<T,FortranStyle> diag() const
    { return SymMatrixView<T,CStyle>::diag(); }

    inline VectorView<T,FortranStyle> diag(int i) const
    { return SymMatrixView<T,CStyle>::diag(i); }

    inline VectorView<T,FortranStyle> diag(int i, int j1, int j2) const
    {
      TMVAssert(j1>0);
      return SymMatrixView<T,CStyle>::diag(i,j1-1,j2); 
    }

    //
    // Modifying Functions
    //

    inline const SymMatrixView<T,FortranStyle>& Zero() const 
    { SymMatrixView<T,CStyle>::Zero(); return *this; }

    inline const SymMatrixView<T,FortranStyle>& SetAllTo(T x) const
    { SymMatrixView<T,CStyle>::SetAllTo(x); return *this; }

    inline const SymMatrixView<T,FortranStyle>& Clip(RealType(T) thresh) const
    { SymMatrixView<T,CStyle>::Clip(thresh); return *this; }

    inline const SymMatrixView<T,FortranStyle>& ConjugateSelf() const
    { SymMatrixView<T,CStyle>::ConjugateSelf(); return *this; }

    inline const SymMatrixView<T,FortranStyle>& TransposeSelf() const
    { SymMatrixView<T,CStyle>::TransposeSelf(); return *this; }

    inline const SymMatrixView<T,FortranStyle>& SetToIdentity(T x=T(1)) const
    { SymMatrixView<T,CStyle>::SetToIdentity(x); return *this; }

    inline const SymMatrixView<T,FortranStyle>& SwapRowsCols(
        int i1, int i2) const
    {
      TMVAssert(i1>0 && i1<=int(this->size()));
      TMVAssert(i2>0 && i2<=int(this->size()));
      SymMatrixView<T,CStyle>::SwapRowsCols(i1-1,i2-1); 
      return *this; 
    }

    inline const SymMatrixView<T,FortranStyle>& PermuteRowsCols(
        const int* p, int i1, int i2) const
    {
      TMVAssert(i1>0 && i1-i2<=0 && i2<=int(this->size()));
      SymMatrixView<T,CStyle>::PermuteRowsCols(p,i1-1,i2); 
      return *this; 
    }

    inline const SymMatrixView<T,FortranStyle>& ReversePermuteRowsCols(
        const int* p, int i1, int i2) const
    {
      TMVAssert(i1>0 && i1-i2<=0 && i2<=int(this->size()));
      SymMatrixView<T,CStyle>::ReversePermuteRowsCols(p,i1-1,i2); 
      return *this; 
    }

    inline const SymMatrixView<T,FortranStyle>& PermuteRowsCols(
        const int* p) const
    { SymMatrixView<T,CStyle>::PermuteRowsCols(p); return *this; }

    inline const SymMatrixView<T,FortranStyle>& ReversePermuteRowsCols(
        const int* p) const
    { SymMatrixView<T,CStyle>::ReversePermuteRowsCols(p); return *this; }

    //
    // SubMatrix
    //

    inline bool OKSubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) const
    {
      return ConstSymMatrixView<T,FortranStyle>(*this).OKSubMatrix(
          i1,i2,j1,j2,istep,jstep); 
    }

    inline bool OKSubVector(
        int i, int j, int istep, int jstep, int n) const
    {
      return ConstSymMatrixView<T,FortranStyle>(*this).OKSubVector(
          i,j,istep,jstep,n); 
    }


    inline bool OKSubSymMatrix(int i1, int i2, int istep) const
    {
      return ConstSymMatrixView<T,FortranStyle>(*this).OKSubSymMatrix(
          i1,i2,istep); 
    }

    inline MatrixView<T,FortranStyle> SubMatrix(
        int i1, int i2, int j1, int j2) const
    {
      TMVAssert(OKSubMatrix(i1,i2,j1,j2,1,1));
      return SymMatrixView<T,CStyle>::SubMatrix(i1-1,i2,j1-1,j2);
    }

    inline MatrixView<T,FortranStyle> SubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) const
    {
      TMVAssert(OKSubMatrix(i1,i2,j1,j2,istep,jstep));
      return SymMatrixView<T,CStyle>::SubMatrix(i1-1,i2-1+istep,
          j1-1,j2-1+jstep,istep,jstep);
    }

    inline VectorView<T,FortranStyle> SubVector(int i, int j,
        int istep, int jstep, int n) const
    {
      TMVAssert(OKSubVector(i,j,istep,jstep,n));
      return SymMatrixView<T,CStyle>::SubVector(i-1,j-1,istep,jstep,n);
    }

    inline SymMatrixView<T,FortranStyle> SubSymMatrix(int i1, int i2) const
    {
      TMVAssert(OKSubSymMatrix(i1,i2,1));
      return SymMatrixView<T,CStyle>::SubSymMatrix(i1-1,i2);
    }

    inline SymMatrixView<T,FortranStyle> SubSymMatrix(int i1, int i2,
        int istep) const
    {
      TMVAssert(OKSubSymMatrix(i1,i2,istep));
      return SymMatrixView<T,CStyle>::SubSymMatrix(i1-1,i2-1+istep,istep);
    }

    inline UpperTriMatrixView<T,FortranStyle> UpperTri(
        DiagType dt = NonUnitDiag) const
    { return SymMatrixView<T,CStyle>::UpperTri(dt); }

    inline LowerTriMatrixView<T,FortranStyle> LowerTri(
        DiagType dt = NonUnitDiag) const
    { return SymMatrixView<T,CStyle>::LowerTri(dt); }

    inline SymMatrixView<T,FortranStyle> View() const
    { return SymMatrixView<T,CStyle>::View(); }

    inline SymMatrixView<T,FortranStyle> Transpose() const
    { return SymMatrixView<T,CStyle>::Transpose(); }

    inline SymMatrixView<T,FortranStyle> Conjugate() const
    { return SymMatrixView<T,CStyle>::Conjugate(); }

    inline SymMatrixView<T,FortranStyle> Adjoint() const
    { return SymMatrixView<T,CStyle>::Adjoint(); }

    inline SymMatrixView<RealType(T),FortranStyle> Real() const
    { return SymMatrixView<T,CStyle>::Real(); }

    inline SymMatrixView<RealType(T),FortranStyle> Imag() const
    { return SymMatrixView<T,CStyle>::Imag(); }

  }; // FortranStyle SymMatrixView

#ifdef XTEST
#ifdef TMVDEBUG
#define XTEST_DEBUG
#endif
#endif

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  class SymMatrix : 
    public GenSymMatrix<T>
  {

  public:

    //
    // Constructors
    //

#define NEW_SIZE(s) itslen((s)*(s)), itsm(new T[itslen]), itss(s)  \
    DEFFIRSTLAST(itsm.get(),itsm.get()+itslen) 

    explicit inline SymMatrix(size_t _size) : NEW_SIZE(_size) 
    { 
      TMVAssert(S==RowMajor || S==ColMajor);
#ifdef TMVDEBUG
      SetAllTo(T(888));
#endif
    }

    inline SymMatrix(size_t _size, const T x) : NEW_SIZE(_size)
    {
      TMVAssert(S==RowMajor || S==ColMajor);
      SetAllTo(x);
    }

    inline SymMatrix(size_t _size, const T* vv) : NEW_SIZE(_size)
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor);
      std::copy(vv,vv+itslen,itsm.get());
    }

    inline SymMatrix(size_t _size, const std::vector<T>& vv) :
      NEW_SIZE(_size)
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor);
      TMVAssert(vv.size() == itslen);
      std::copy(vv.begin(),vv.end(),itsm.get());
    }

    inline SymMatrix(const SymMatrix<T,U,S,I>& rhs) : NEW_SIZE(rhs.size())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor);
      std::copy(rhs.cptr(),rhs.cptr()+itslen,itsm.get());
    }

    template <IndexStyle I2> 
    inline SymMatrix(const SymMatrix<T,U,S,I2>& rhs) : NEW_SIZE(rhs.size())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor);
      std::copy(rhs.cptr(),rhs.cptr()+itslen,itsm.get());
    }

    inline SymMatrix(const GenSymMatrix<RealType(T)>& rhs) :
      NEW_SIZE(rhs.size())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor);
      if (rhs.issym()) 
        rhs.AssignToS(View());
      else {
        if (U == Upper)
          UpperTri() = rhs.UpperTri();
        else
          LowerTri() = rhs.LowerTri();
      }
    }

    inline SymMatrix(const GenSymMatrix<ComplexType(T)>& rhs) :
      NEW_SIZE(rhs.size())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(IsComplex(T()));
      TMVAssert(S==RowMajor || S==ColMajor);
      if (rhs.issym()) 
        rhs.AssignToS(View());
      else {
        if (U == Upper)
          UpperTri() = rhs.UpperTri();
        else
          LowerTri() = rhs.LowerTri();
      }
    }

    template <class T2> 
    inline SymMatrix(const GenSymMatrix<T2>& rhs) :
      NEW_SIZE(rhs.size())
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor);
      if (U == Upper)
        UpperTri() = rhs.UpperTri();
      else
        LowerTri() = rhs.LowerTri();
    }

    template <IndexStyle I2> 
    inline explicit SymMatrix(const Matrix<T,S,I2>& rhs) :
      NEW_SIZE(rhs.rowsize())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor);
      TMVAssert(rhs.IsSquare());
      std::copy(rhs.cptr(),rhs.cptr()+itslen,itsm.get());
    }

    template <class T2> 
    inline SymMatrix(const GenMatrix<T2>& rhs) :
      NEW_SIZE(rhs.rowsize())
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor);
      TMVAssert(rhs.IsSquare());
      if (U == Upper)
        UpperTri() = rhs.UpperTri();
      else
        LowerTri() = rhs.LowerTri();
    }

    inline SymMatrix(const AssignableToSymMatrix<RealType(T)>& m2) :
      NEW_SIZE(m2.size())
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor);
      TMVAssert(m2.issym());
      m2.AssignToS(View());
    }

    inline SymMatrix(const AssignableToSymMatrix<ComplexType(T)>& m2) :
      NEW_SIZE(m2.size())
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(IsComplex(T()));
      TMVAssert(S==RowMajor || S==ColMajor);
      TMVAssert(m2.issym());
      m2.AssignToS(View());
    }

    inline explicit SymMatrix(const GenDiagMatrix<RealType(T)>& m2) :
      NEW_SIZE(m2.size())
    { 
      TMVAssert(size() == m2.size());
      Zero();
      m2.AssignToD(DiagMatrixViewOf(diag()));
    }

    inline explicit SymMatrix(const GenDiagMatrix<ComplexType(T)>& m2) :
      NEW_SIZE(m2.size())
    { 
      TMVAssert(IsComplex(T()));
      TMVAssert(size() == m2.size());
      Zero();
      m2.AssignToD(DiagMatrixViewOf(diag()));
    }

    template <class T2> 
    inline SymMatrix(const OProdVV<T,T2,T2>& opvv) :
      NEW_SIZE(opvv.colsize())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(opvv.colsize() == opvv.rowsize());
      TMVAssert(opvv.GetV1().SameAs(opvv.GetV2().View()));
      Rank1Update<false>(opvv.GetX(),opvv.GetV1(),View());
    }

    template <class T2> 
    inline SymMatrix(const ProdMM<T,T2,T2>& pmm) :
      NEW_SIZE(pmm.colsize())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(pmm.colsize() == pmm.rowsize());
      TMVAssert(pmm.GetM1().SameAs(pmm.GetM2().Transpose()));
      RankKUpdate<false>(pmm.GetX(),pmm.GetM1(),View());
    }

    template <class T2> 
    inline SymMatrix(const ProdUL<T,T2,T2>& pmm) :
      NEW_SIZE(pmm.colsize())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(pmm.colsize() == pmm.rowsize());
      TMVAssert(pmm.GetM1().SameAs(pmm.GetM2().Transpose()));
      RankKUpdate<false>(pmm.GetX(),pmm.GetM1(),View());
    }

    template <class T2> 
    inline SymMatrix(const ProdLU<T,T2,T2>& pmm) :
      NEW_SIZE(pmm.colsize())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(pmm.colsize() == pmm.rowsize());
      TMVAssert(pmm.GetM1().SameAs(pmm.GetM2().Transpose()));
      RankKUpdate<false>(pmm.GetX(),pmm.GetM1(),View());
    }

#undef NEW_SIZE

    virtual inline ~SymMatrix() 
    {
#ifdef TMVDEBUG
      SetAllTo(T(999));
#endif
    }

    //
    // Op=
    //

    inline SymMatrix<T,U,S,I>& operator=(const SymMatrix<T,U,S,I>& m2)
    { 
      TMVAssert(size() == m2.size());
      if (&m2 != this) 
        std::copy(m2.cptr(),m2.cptr()+itslen,itsm.get());
      return *this;
    }

    template <IndexStyle I2>
    inline SymMatrix<T,U,S,I>& operator=(const SymMatrix<T,U,S,I2>& m2)
    { 
      TMVAssert(size() == m2.size());
      if (&m2 != this) 
        std::copy(m2.cptr(),m2.cptr()+itslen,itsm.get());
      return *this;
    }

    inline SymMatrix<T,U,S,I>& operator=(const GenSymMatrix<RealType(T)>& m2)
    { 
      TMVAssert(size() == m2.size());
      TMVAssert(m2.issym());
      m2.AssignToS(View());
      return *this; 
    }

    inline SymMatrix<T,U,S,I>& operator=(
        const GenSymMatrix<ComplexType(T)>& m2)
    { 
      TMVAssert(IsComplex(T()));
      TMVAssert(size() == m2.size());
      TMVAssert(m2.issym());
      m2.AssignToS(View());
      return *this; 
    }

    template <class T2> 
    inline SymMatrix<T,U,S,I>& operator=(const GenSymMatrix<T2>& m2)
    { 
      TMVAssert(size() == m2.size());
      TMVAssert(IsReal(T2()) || m2.issym());
      UpperTri() = m2.UpperTri(); 
      return *this; 
    }

    inline SymMatrix<T,U,S,I>& operator=(T x) 
    { return SetToIdentity(x); }

    inline SymMatrix<T,U,S,I>& operator=(
        const AssignableToSymMatrix<RealType(T)>& m2) 
    { 
      TMVAssert(size() == m2.size());
      TMVAssert(m2.issym());
      m2.AssignToS(View());
      return *this;
    }

    inline SymMatrix<T,U,S,I>& operator=(
        const AssignableToSymMatrix<ComplexType(T)>& m2) 
    { 
      TMVAssert(IsComplex(T()));
      TMVAssert(size() == m2.size());
      TMVAssert(m2.issym());
      m2.AssignToS(View());
      return *this;
    }

    inline SymMatrix<T,U,S,I>& operator=(
        const GenDiagMatrix<RealType(T)>& m2) 
    { 
      TMVAssert(size() == m2.size());
      View() = m2;
      return *this;
    }

    inline SymMatrix<T,U,S,I>& operator=(
        const GenDiagMatrix<ComplexType(T)>& m2) 
    { 
      TMVAssert(IsComplex(T()));
      TMVAssert(size() == m2.size());
      View() = m2;
      return *this;
    }

    template <class T2> 
    inline SymMatrix<T,U,S,I>& operator=(const OProdVV<T,T2,T2>& opvv)
    {
      TMVAssert(size() == opvv.colsize());
      TMVAssert(size() == opvv.rowsize());
      TMVAssert(opvv.GetV1().SameAs(opvv.GetV2().View()));
      Rank1Update<false>(opvv.GetX(),opvv.GetV1(),View());
      return *this;
    }

    template <class T2> 
    inline SymMatrix<T,U,S,I>& operator=(const ProdMM<T,T2,T2>& pmm) 
    {
      TMVAssert(size() == pmm.colsize());
      TMVAssert(size() == pmm.rowsize());
      TMVAssert(pmm.GetM1().SameAs(pmm.GetM2().Transpose()));
      RankKUpdate<false>(pmm.GetX(),pmm.GetM1(),View());
      return *this;
    }

    template <class T2> 
    inline SymMatrix<T,U,S,I>& operator=(const ProdUL<T,T2,T2>& pmm)
    {
      TMVAssert(size() == pmm.colsize());
      TMVAssert(size() == pmm.rowsize());
      TMVAssert(pmm.GetM1().SameAs(pmm.GetM2().Transpose()));
      RankKUpdate<false>(pmm.GetX(),pmm.GetM1(),View());
      return *this;
    }

    template <class T2> 
    inline SymMatrix<T,U,S,I>& operator=(const ProdLU<T,T2,T2>& pmm) 
    {
      TMVAssert(size() == pmm.colsize());
      TMVAssert(size() == pmm.rowsize());
      TMVAssert(pmm.GetM1().SameAs(pmm.GetM2().Transpose()));
      RankKUpdate<false>(pmm.GetX(),pmm.GetM1(),View());
      return *this;
    }


    //
    // Access
    //

    typedef T& reference;

    inline T operator()(int i, int j) const
    { 
      if (I==CStyle) {
        TMVAssert(i>=0 && i<int(size()));
        TMVAssert(j>=0 && j<int(size()));
        return cref(i,j); 
      } else {
        TMVAssert(i>0 && i<=int(size()));
        TMVAssert(j>0 && j<=int(size()));
        return cref(i-1,j-1); 
      }
    }

    inline T& operator()(int i, int j) 
    { 
      if (I==CStyle) {
        TMVAssert(i>=0 && i<int(size()));
        TMVAssert(j>=0 && j<int(size()));
        return ref(i,j); 
      } else {
        TMVAssert(i>0 && i<=int(size()));
        TMVAssert(j>0 && j<=int(size()));
        return ref(i-1,j-1); 
      }
    }

    inline ConstVectorView<T,I> row(int i, int j1, int j2) const 
    { 
      if (I==FortranStyle) {
        TMVAssert(i>0 && i<=int(size())); --i;
        TMVAssert(j1>0 && j1-j2<=0 && j2<=int(size())); --j1;
      } else {
        TMVAssert(i>=0 && i<int(size()));
        TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size()));
      }
      TMVAssert( i-j1<=0 || j2-i<=1 );
      if ((U==Upper && i-j1<=0) || (U==Lower && j2-i<=1))
        return ConstVectorView<T,I>(itsm.get()+i*stepi()+j1*stepj(),
            j2-j1,stepj(),NonConj); 
      else
        return ConstVectorView<T,I>(itsm.get()+i*stepj()+j1*stepi(),
            j2-j1,stepi(),NonConj); 
    }

    inline ConstVectorView<T,I> col(int j, int i1, int i2) const
    {
      if (I==FortranStyle) { 
        TMVAssert(j>0 && j<=int(size())); --j;
        TMVAssert(i1>0 && i1-i2<=0 && i2<=int(size())); --i1;
      } else {
        TMVAssert(j>=0 && j<int(size()));
        TMVAssert(i1>=0 && i1-i2<=0 && i2<=int(size()));
      }
      TMVAssert( j-i1<=0 || i2-j<=1 );
      if ((U==Upper && i2-j<=1) || (U==Lower && j-i1<=0))
        return ConstVectorView<T,I>(itsm.get()+i1*stepi()+j*stepj(),
            i2-i1,stepi(),NonConj); 
      else 
        return ConstVectorView<T,I>(itsm.get()+i1*stepj()+j*stepi(),
            i2-i1,stepj(),NonConj); 
    }

    inline ConstVectorView<T,I> diag() const
    { return ConstVectorView<T,I>(itsm.get(),size(),stepi()+stepj(),NonConj); }

    inline ConstVectorView<T,I> diag(int i) const
    {
      TMVAssert(i>=-int(size()) && i<=int(size())); 
      i = std::abs(i);
      return ConstVectorView<T,I>(itsm.get()+i*(U==Upper?stepj():stepi()),
          size()-i,stepi()+stepj(),NonConj); 
    }

    inline ConstVectorView<T,I> diag(int i, int j1, int j2) const
    {
      TMVAssert(i>=-int(size()) && i<=int(size())); 
      if (I==FortranStyle) {
        TMVAssert(j1>0 && j1-j2<=0 && j2<=int(size())-std::abs(i)); --j1;
      } else {
        TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size())-std::abs(i));
      }
      i = std::abs(i);
      const int ds = stepi()+stepj();
      return ConstVectorView<T,I>(itsm.get()+i*(U==Upper?stepj():stepi())+j1*ds,
          j2-j1,ds,NonConj);
    }

    inline VectorView<T,I> row(int i, int j1, int j2)
    { 
      if (I==FortranStyle) {
        TMVAssert(i>0 && i<=int(size())); --i;
        TMVAssert(j1>0 && j1-j2<=0 && j2<=int(size())); --j1;
      } else {
        TMVAssert(i>=0 && i<int(size()));
        TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size()));
      }
      TMVAssert( i-j1<=0 || j2-i<=1 );
      if ((U==Upper && i-j1<=0) || (U==Lower && j2-i<=1))
        return VectorView<T,I>(itsm.get()+i*stepi()+j1*stepj(),
            j2-j1,stepj(),NonConj FIRSTLAST); 
      else
        return VectorView<T,I>(itsm.get()+i*stepj()+j1*stepi(),
            j2-j1,stepi(),NonConj FIRSTLAST); 
    }

    inline VectorView<T,I> col(int j, int i1, int i2)
    {
      if (I==FortranStyle) { 
        TMVAssert(j>0 && j<=int(size())); --j;
        TMVAssert(i1>0 && i1-i2<=0 && i2<=int(size())); --i1;
      } else {
        TMVAssert(j>=0 && j<int(size()));
        TMVAssert(i1>=0 && i1-i2<=0 && i2<=int(size()));
      }
      TMVAssert( j-i1<=0 || i2-j<=1 );
      if ((U==Upper && i2-j<=1) || (U==Lower && j-i1<=0))
        return VectorView<T,I>(itsm.get()+i1*stepi()+j*stepj(),
            i2-i1,stepi(),NonConj FIRSTLAST); 
      else 
        return VectorView<T,I>(itsm.get()+i1*stepj()+j*stepi(),
            i2-i1,stepj(),NonConj FIRSTLAST); 
    }

    inline VectorView<T,I> diag()
    {
      return VectorView<T,I>(itsm.get(),size(),stepi()+stepj(),NonConj 
          FIRSTLAST);
    }

    inline VectorView<T,I> diag(int i) 
    {
      TMVAssert(i>=-int(size()) && i<=int(size())); 
      i = std::abs(i);
      TMVAssert(i<=int(size())); 
      return VectorView<T,I>(itsm.get()+i*(U==Upper?stepj():stepi()),
          size()-i,stepi()+stepj(),NonConj FIRSTLAST); 
    }

    inline VectorView<T,I> diag(int i, int j1, int j2) 
    {
      TMVAssert(i>=-int(size()) && i<=int(size())); 
      if (I==FortranStyle) {
        TMVAssert(j1>0 && j1-j2<=0 && j2<=int(size())-std::abs(i)); --j1;
      } else {
        TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size())-std::abs(i));
      }
      i = std::abs(i);
      const int ds = stepi()+stepj();
      return VectorView<T,I>(itsm.get()+i*(U==Upper?stepj():stepi())+j1*ds,
          j2-j1,ds,NonConj FIRSTLAST);
    }


    //
    // Modifying Functions
    //

    inline SymMatrix<T,U,S,I>& Zero() 
    {
      std::fill_n(itsm.get(),itslen,T(0));
      return *this;
    }

    inline SymMatrix<T,U,S,I>& SetAllTo(T x) 
    { UpperTri().SetAllTo(x); return *this; }

    inline SymMatrix<T,U,S,I>& Clip(RealType(T) thresh) 
    { UpperTri().Clip(thresh); return *this; }

    inline SymMatrix<T,U,S,I>& ConjugateSelf() 
    { if (IsComplex(T())) UpperTri().ConjugateSelf(); return *this; }

    inline SymMatrix<T,U,S,I>& TransposeSelf() 
    { return *this; }

    inline SymMatrix<T,U,S,I>& SetToIdentity(T x=T(1)) 
    { Zero(); diag().SetAllTo(x); return *this; }

    inline SymMatrix<T,U,S,I>& SwapRowsCols(int i1, int i2)
    { View().SwapRowsCols(i1,i2); return *this; }

    inline SymMatrix<T,U,S,I>& PermuteRowsCols(const int* p,
        int i1, int i2)
    { View().PermuteRowsCols(p,i1,i2); return *this; }

    inline SymMatrix<T,U,S,I>& ReversePermuteRowsCols(const int* p,
        int i1, int i2)
    { View().ReversePermuteRowsCols(p,i1,i2); return *this; }

    inline SymMatrix<T,U,S,I>& PermuteRowsCols(const int* p)
    { View().PermuteRowsCols(p); return *this; }

    inline SymMatrix<T,U,S,I>& ReversePermuteRowsCols(const int* p)
    { View().ReversePermuteRowsCols(p); return *this; }


    //
    // SubMatrix
    //

    inline ConstMatrixView<T,I> SubMatrix(
        int i1, int i2, int j1, int j2) const
    {
      TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,1,1));
      if (I==FortranStyle) { --i1; --j1; }
      TMVAssert(i2-1<=j1 || j2-1<=i1);
      if ((U==Upper && i2-1<=j1) || (U==Lower && j2-1<=i1))
        return ConstMatrixView<T,I>(itsm.get()+i1*stepi()+j1*stepj(),
            i2-i1,j2-j1,stepi(),stepj(),S,NonConj);
      else
        return ConstMatrixView<T,I>(itsm.get()+i1*stepj()+j1*stepi(),
            i2-i1,j2-j1,stepj(),stepi(),TransOf(S),NonConj);
    }

    inline ConstMatrixView<T,I> SubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) const
    {
      TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,istep,jstep));
      if (I==FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
      const StorageType newstor = S==RowMajor ?
      jstep == 1 ? RowMajor : NoMajor :
      istep == 1 ? ColMajor : NoMajor;
      TMVAssert(i2-istep<=j1 || j2-jstep<=i1);
      if ((U==Upper && i2-istep<=j1) || (U==Lower && j2-jstep<=i1))
        return ConstMatrixView<T,I>(itsm.get()+i1*stepi()+j1*stepj(),
            (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
            newstor,NonConj);
      else
        return ConstMatrixView<T,I>(itsm.get()+i1*stepj()+j1*stepi(),
            (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
            TransOf(newstor),NonConj);
    }

    inline ConstVectorView<T,I> SubVector(int i, int j,
        int istep, int jstep, int n) const
    {
      TMVAssert(View().OKSubVector(i,j,istep,jstep,n));
      if (I==FortranStyle) { --i; --j; }
      if ((U==Upper && i-j<=0) || (U==Lower && j-i<=0))
        return ConstVectorView<T,I>(itsm.get()+i*stepi()+j*stepj(),n,
            istep*stepi()+jstep*stepj(),NonConj);
      else
        return ConstVectorView<T,I>(itsm.get()+i*stepj()+j*stepi(),n,
            istep*stepj()+jstep*stepi(),NonConj);
    }

    inline ConstSymMatrixView<T,I> SubSymMatrix(int i1, int i2) const
    {
      TMVAssert(View().OKSubSymMatrix(i1,i2,1));
      if (I==FortranStyle) { --i1; }
      return ConstSymMatrixView<T,I>(itsm.get()+i1*(stepi()+stepj()),i2-i1,
          stepi(),stepj(),Sym,U,S,NonConj);
    }

    inline ConstSymMatrixView<T,I> SubSymMatrix(int i1, int i2,
        int istep) const
    {
      TMVAssert(View().OKSubSymMatrix(i1,i2,istep));
      if (I==FortranStyle) { --i1; i2+=istep-1; }
      return ConstSymMatrixView<T,I>(itsm.get()+i1*(stepi()+stepj()),
          (i2-i1)/istep,istep*stepi(),istep*stepj(),Sym,U,
          istep==1 ? S : NoMajor, NonConj);
    }

    inline ConstUpperTriMatrixView<T,I> UpperTri(
        DiagType dt = NonUnitDiag) const
    {
      return U==Upper ? 
      ConstUpperTriMatrixView<T,I>(itsm.get(),size(),
          stepi(),stepj(),dt,S,NonConj) :
      ConstUpperTriMatrixView<T,I>(itsm.get(),size(),
          stepj(),stepi(),dt,TransOf(S),NonConj);
    }

    inline ConstLowerTriMatrixView<T,I> LowerTri(
        DiagType dt = NonUnitDiag) const
    {
      return U==Lower ? 
      ConstLowerTriMatrixView<T,I>(itsm.get(),size(),
          stepi(),stepj(),dt,S,NonConj) :
      ConstLowerTriMatrixView<T,I>(itsm.get(),size(),
          stepj(),stepi(),dt,TransOf(S),NonConj);
    }

    inline MatrixView<T,I> SubMatrix(int i1, int i2, int j1, int j2)
    {
      TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,1,1));
      if (I==FortranStyle) { --i1; --j1; }
      TMVAssert(i2-1<=j1 || j2-1<=i1);
      if ((U==Upper && i2-1<=j1) || (U==Lower && j2-1<=i1))
        return MatrixView<T,I>(itsm.get()+i1*stepi()+j1*stepj(),
            i2-i1,j2-j1,stepi(),stepj(),S,NonConj FIRSTLAST);
      else
        return MatrixView<T,I>(itsm.get()+i1*stepj()+j1*stepi(),
            i2-i1,j2-j1,stepj(),stepi(),TransOf(S),NonConj FIRSTLAST);
    }

    inline MatrixView<T,I> SubMatrix(int i1, int i2, int j1, int j2,
        int istep, int jstep)
    {
      const StorageType newstor = S==RowMajor ?
      jstep == 1 ? RowMajor : NoMajor :
      istep == 1 ? ColMajor : NoMajor;
      TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,istep,jstep));
      if (I==FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
      TMVAssert(i2-istep<=j1 || j2-jstep<=i1);
      if ((U==Upper && i2-istep<=j1) || (U==Lower && j2-jstep<=i1))
        return MatrixView<T,I>(itsm.get()+i1*stepi()+j1*stepj(),
            (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
            newstor,NonConj FIRSTLAST);
      else
        return MatrixView<T,I>(itsm.get()+i1*stepj()+j1*stepi(),
            (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
            TransOf(newstor),NonConj FIRSTLAST);
    }

    inline VectorView<T,I> SubVector(int i, int j,
        int istep, int jstep, int n)
    {
      TMVAssert(View().OKSubVector(i,j,istep,jstep,n));
      if (I==FortranStyle) { --i; --j; }
      if ((U==Upper && i-j<=0) || (U==Lower && j-i<=0))
        return VectorView<T,I>(itsm.get()+i*stepi()+j*stepj(),n,
            istep*stepi()+jstep*stepj(),NonConj FIRSTLAST);
      else
        return VectorView<T,I>(itsm.get()+i*stepj()+j*stepi(),n,
            istep*stepj()+jstep*stepi(),NonConj FIRSTLAST);
    }

    inline SymMatrixView<T,I> SubSymMatrix(int i1, int i2)
    {
      TMVAssert(View().OKSubSymMatrix(i1,i2,1));
      if (I==FortranStyle) { --i1; }
      return SymMatrixView<T,I>(itsm.get()+i1*(stepi()+stepj()),i2-i1,
          stepi(),stepj(),Sym,U,S,NonConj FIRSTLAST);
    }

    inline SymMatrixView<T,I> SubSymMatrix(int i1, int i2, int istep) 
    {
      TMVAssert(View().OKSubSymMatrix(i1,i2,istep));
      if (I==FortranStyle) { --i1; i2+=istep-1; }
      return SymMatrixView<T,I>(itsm.get()+i1*(stepi()+stepj()),(i2-i1)/istep,
          istep*stepi(),istep*stepj(),Sym,U,
          istep==1 ? S : NoMajor,NonConj FIRSTLAST);
    }

    inline UpperTriMatrixView<T,I> UpperTri(DiagType dt = NonUnitDiag)
    {
      return U==Upper ? 
      UpperTriMatrixView<T,I>(itsm.get(),size(),
          stepi(),stepj(),dt,S,NonConj FIRSTLAST) :
      UpperTriMatrixView<T,I>(itsm.get(),size(),
          stepj(),stepi(),dt,TransOf(S),NonConj FIRSTLAST);
    }

    inline LowerTriMatrixView<T,I> LowerTri(DiagType dt = NonUnitDiag)
    {
      return U==Lower ? 
      LowerTriMatrixView<T,I>(itsm.get(),size(),
          stepi(),stepj(),dt,S,NonConj FIRSTLAST) :
      LowerTriMatrixView<T,I>(itsm.get(),size(),
          stepj(),stepi(),dt,TransOf(S),NonConj FIRSTLAST);
    }

    inline ConstSymMatrixView<T,I> View() const
    { 
      return ConstSymMatrixView<T,I>(itsm.get(),size(),stepi(),stepj(),
          Sym,U,S,NonConj);
    }

    inline ConstSymMatrixView<T,I> Transpose() const
    {
      return ConstSymMatrixView<T,I>(itsm.get(),size(),stepj(),stepi(),
          Sym,UTransOf(U),TransOf(S),NonConj);
    }

    inline ConstSymMatrixView<T,I> Conjugate() const
    { 
      return ConstSymMatrixView<T,I>(itsm.get(),size(),stepi(),stepj(),
          Sym,U,S,ConjOf(T,NonConj));
    }

    inline ConstSymMatrixView<T,I> Adjoint() const
    {
      return ConstSymMatrixView<T,I>(itsm.get(),size(),stepj(),stepi(),
          Sym,UTransOf(U),TransOf(S),ConjOf(T,NonConj));
    }

    inline ConstSymMatrixView<RealType(T),I> Real() const
    {
      return ConstSymMatrixView<RealType(T)>(
          reinterpret_cast<const RealType(T)*>(itsm.get()),size(),
          IsReal(T()) ? stepi() : 2*stepi(),
          IsReal(T()) ? stepj() : 2*stepj(),
          Sym,U,IsReal(T())?S:NoMajor,NonConj);
    }

    inline ConstSymMatrixView<RealType(T),I> Imag() const
    {
      TMVAssert(IsComplex(T()));
      return ConstSymMatrixView<RealType(T)>(
          reinterpret_cast<const RealType(T)*>(itsm.get())+1,size(),
          2*stepi(),2*stepj(),Sym,U,NoMajor,NonConj);
    } 

    inline SymMatrixView<T,I> View() 
    { 
      return SymMatrixView<T,I>(itsm.get(),size(),stepi(),stepj(),
          Sym,U,S,NonConj FIRSTLAST);
    }

    inline SymMatrixView<T,I> Transpose() 
    {
      return SymMatrixView<T,I>(itsm.get(),size(),stepj(),stepi(),
          Sym,UTransOf(U),TransOf(S),NonConj FIRSTLAST);
    }

    inline SymMatrixView<T,I> Conjugate() 
    { 
      return SymMatrixView<T,I>(itsm.get(),size(),stepi(),stepj(),
          Sym,U,S,ConjOf(T,NonConj) FIRSTLAST);
    }

    inline SymMatrixView<T,I> Adjoint() 
    {
      return SymMatrixView<T,I>(itsm.get(),size(),stepj(),stepi(),
          Sym,UTransOf(U),TransOf(S),ConjOf(T,NonConj) FIRSTLAST);
    }

    inline SymMatrixView<RealType(T),I> Real()
    {
      return SymMatrixView<RealType(T)>(
          reinterpret_cast<RealType(T)*>(itsm.get()),size(),
          IsReal(T()) ? stepi() : 2*stepi(),
          IsReal(T()) ? stepj() : 2*stepj(),
          Sym,U,IsReal(T())?S:NoMajor,NonConj
#ifdef TMVFLDEBUG
          ,reinterpret_cast<const RealType(T)*>(first)
          ,reinterpret_cast<const RealType(T)*>(last)
#endif
          );
    }

    inline SymMatrixView<RealType(T),I> Imag()
    {
      TMVAssert(IsComplex(T()));
      return SymMatrixView<RealType(T)>(
          reinterpret_cast<RealType(T)*>(itsm.get())+1,size(),
          2*stepi(),2*stepj(),Sym,U,NoMajor,NonConj
#ifdef TMVFLDEBUG
          ,reinterpret_cast<const RealType(T)*>(first)+1
          ,reinterpret_cast<const RealType(T)*>(last)+1
#endif
          );
    } 

    inline size_t size() const { return itss; }
    inline const T* cptr() const { return itsm.get(); }
    inline T* ptr() { return itsm.get(); }
    inline int stepi() const { return S==RowMajor ? itss : 1; }
    inline int stepj() const { return S==RowMajor ? 1 : itss; }
    inline SymType sym() const { return Sym; }
    inline UpLoType uplo() const { return U; }
    inline StorageType stor() const { return S; }
    inline ConjType ct() const { return NonConj; }
    inline bool isrm() const { return S==RowMajor; }
    inline bool iscm() const { return S==ColMajor; }
    inline bool isconj() const { return false; }
    inline bool isherm() const { return IsReal(T()); }
    inline bool issym() const { return true; }
    inline bool isupper() const { return U == Upper; }

    inline T& ref(int i, int j)
    {
      if ((U==Upper && i <= j) || (U==Lower && i>=j)) 
        if (S == RowMajor) return itsm.get()[i*itss + j];
        else return itsm.get()[j*itss + i];
      else 
        if (S == RowMajor) return itsm.get()[j*itss + i];
        else return itsm.get()[i*itss + j];
    }

    inline T cref(int i, int j) const 
    {
      if ((U==Upper && i <= j) || (U==Lower && i>=j)) 
        if (S == RowMajor) return itsm.get()[i*itss + j];
        else return itsm.get()[j*itss + i];
      else 
        if (S == RowMajor) return itsm.get()[j*itss + i];
        else return itsm.get()[i*itss + j];
    }

  protected :

    const size_t itslen;
    auto_array<T> itsm;
    const size_t itss;

#ifdef TMVFLDEBUG
  public:
    const T*const first;
    const T*const last;
#endif

  }; // SymMatrix

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  class HermMatrix : 
    public GenSymMatrix<T>
  {

  public:

    //
    // Constructors
    //

#define NEW_SIZE(s) itslen((s)*(s)), itsm(new T[itslen]), itss(s)  \
    DEFFIRSTLAST(itsm.get(),itsm.get()+itslen)

    explicit inline HermMatrix(size_t _size) : NEW_SIZE(_size) 
    { 
      TMVAssert(S==RowMajor || S==ColMajor);
#ifdef TMVDEBUG
      SetAllTo(RealType(T)(888));
#endif
    }

    HermMatrix(size_t _size, const RealType(T) x) : NEW_SIZE(_size)
    {
      TMVAssert(S==RowMajor || S==ColMajor);
      SetAllTo(x);
    }

    HermMatrix(size_t _size, const T* vv) : NEW_SIZE(_size)
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor);
      std::copy(vv,vv+itslen,itsm.get());
      TMVAssert(this->HermOK());
    }

    HermMatrix(size_t _size, const std::vector<T>& vv) : NEW_SIZE(_size)
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor);
      TMVAssert(vv.size() == itslen);
      std::copy(vv.begin(),vv.end(),itsm.get());
      TMVAssert(this->HermOK());
    }

    HermMatrix(const HermMatrix<T,U,S,I>& rhs) : NEW_SIZE(rhs.size())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor);
      std::copy(rhs.cptr(),rhs.cptr()+itslen,itsm.get());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    template <IndexStyle I2> 
    HermMatrix(const HermMatrix<T,U,S,I2>& rhs) :
      NEW_SIZE(rhs.size())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor);
      std::copy(rhs.cptr(),rhs.cptr()+itslen,itsm.get());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    HermMatrix(const GenSymMatrix<RealType(T)>& rhs) : NEW_SIZE(rhs.size())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor);
      if (rhs.isherm()) rhs.AssignToS(View());
      else if (U == Upper) UpperTri() = rhs.UpperTri();
      else LowerTri() = rhs.LowerTri();
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    HermMatrix(const GenSymMatrix<ComplexType(T)>& rhs) : 
      NEW_SIZE(rhs.size())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(IsComplex(T()));
      TMVAssert(S==RowMajor || S==ColMajor);
      if (rhs.isherm())
        rhs.AssignToS(View());
      else {
        if (U == Upper) UpperTri() = rhs.UpperTri();
        else LowerTri() = rhs.LowerTri();
        diag().Imag().Zero();
      }
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    template <class T2> 
    inline HermMatrix(const GenSymMatrix<T2>& rhs) :
      NEW_SIZE(rhs.size())
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor);
      if (U == Upper) UpperTri() = rhs.UpperTri();
      else LowerTri() = rhs.LowerTri();
      if (IsComplex(T()) && IsComplex(T2()) && rhs.issym()) 
        diag().Imag().Zero();
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    template <IndexStyle I2> 
    inline explicit HermMatrix(const Matrix<T,S,I2>& rhs) : 
      NEW_SIZE(rhs.rowsize())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor);
      TMVAssert(rhs.IsSquare());
      std::copy(rhs.cptr(),rhs.cptr()+itslen,itsm.get());
      if (IsComplex(T())) diag().Imag().Zero();
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    template <class T2> 
    inline HermMatrix(const GenMatrix<T2>& rhs) :
      NEW_SIZE(rhs.rowsize())
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor);
      TMVAssert(rhs.IsSquare());
      if (U == Upper)
        UpperTri() = rhs.UpperTri();
      else
        LowerTri() = rhs.LowerTri();
      if (IsComplex(T()) && IsComplex(T2())) diag().Imag().Zero();
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    inline HermMatrix(const AssignableToSymMatrix<RealType(T)>& m2) :
      NEW_SIZE(m2.size())
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor);
      TMVAssert(m2.isherm());
      m2.AssignToS(View());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    inline HermMatrix(const AssignableToSymMatrix<ComplexType(T)>& m2) :
      NEW_SIZE(m2.size())
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(IsComplex(T()));
      TMVAssert(S==RowMajor || S==ColMajor);
      TMVAssert(m2.isherm());
      m2.AssignToS(View());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    inline explicit HermMatrix(const GenDiagMatrix<RealType(T)>& m2) :
      NEW_SIZE(m2.size())
    { 
      TMVAssert(S==RowMajor || S==ColMajor);
      Zero();
      m2.AssignToD(DiagMatrixViewOf(diag()));
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    inline explicit HermMatrix(const GenDiagMatrix<ComplexType(T)>& m2) :
      NEW_SIZE(m2.size())
    { 
      TMVAssert(S==RowMajor || S==ColMajor);
      Zero();
      m2.AssignToD(DiagMatrixViewOf(diag().Real()));
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    template <class T2> 
    inline HermMatrix(const OProdVV<T,T2,T2>& opvv) :
      NEW_SIZE(opvv.colsize())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(opvv.colsize() == opvv.rowsize());
      TMVAssert(opvv.GetV1().SameAs(opvv.GetV2().Conjugate()));
      TMVAssert(IMAG(opvv.GetX()) == RealType(T)(0));
      Rank1Update<false>(opvv.GetX(),opvv.GetV1(),View());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    template <class T2> 
    inline HermMatrix(const ProdMM<T,T2,T2>& pmm) :
      NEW_SIZE(pmm.colsize())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(pmm.colsize() == pmm.rowsize());
      TMVAssert(pmm.GetM1().SameAs(pmm.GetM2().Adjoint()));
      TMVAssert(IMAG(pmm.GetX()) == RealType(T)(0));
      RankKUpdate<false>(pmm.GetX(),pmm.GetM1(),View());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    template <class T2> 
    inline HermMatrix(const ProdUL<T,T2,T2>& pmm) :
      NEW_SIZE(pmm.colsize())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(pmm.colsize() == pmm.rowsize());
      TMVAssert(pmm.GetM1().SameAs(pmm.GetM2().Adjoint()));
      TMVAssert(IMAG(pmm.GetX()) == RealType(T)(0));
      RankKUpdate<false>(pmm.GetX(),pmm.GetM1(),View());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    template <class T2> 
    inline HermMatrix(const ProdLU<T,T2,T2>& pmm) :
      NEW_SIZE(pmm.colsize())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(pmm.colsize() == pmm.rowsize());
      TMVAssert(pmm.GetM1().SameAs(pmm.GetM2().Adjoint()));
      TMVAssert(IMAG(pmm.GetX()) == RealType(T)(0));
      RankKUpdate<false>(pmm.GetX(),pmm.GetM1(),View());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }


#undef NEW_SIZE

    virtual inline ~HermMatrix() 
    {
#ifdef TMVDEBUG
      SetAllTo(T(999));
#endif
    }

    //
    // Op=
    //

    inline HermMatrix<T,U,S,I>& operator=(const HermMatrix<T,U,S,I>& m2)
    { 
      TMVAssert(m2.size() == size());
      if (&m2 != this) 
        std::copy(m2.cptr(),m2.cptr()+itslen,itsm.get());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    template <IndexStyle I2>
    inline HermMatrix<T,U,S,I>& operator=(const HermMatrix<T,U,S,I2>& m2)
    { 
      TMVAssert(m2.size() == size());
      if (&m2 != this) 
        std::copy(m2.cptr(),m2.cptr()+itslen,itsm.get());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    inline HermMatrix<T,U,S,I>& operator=(const GenSymMatrix<RealType(T)>& m2)
    {
      TMVAssert(m2.size() == size());
      TMVAssert(m2.isherm());
      m2.AssignToS(View());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this; 
    }

    inline HermMatrix<T,U,S,I>& operator=(
        const GenSymMatrix<ComplexType(T)>& m2)
    {
      TMVAssert(IsComplex(T()));
      TMVAssert(m2.size() == size());
      TMVAssert(m2.isherm());
      m2.AssignToS(View());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this; 
    }

    template <class T2> 
    inline HermMatrix<T,U,S,I>& operator=(const GenSymMatrix<T2>& m2)
    {
      TMVAssert(m2.size() == size());
      TMVAssert(m2.isherm());
      UpperTri() = m2.UpperTri(); 
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this; 
    }

    inline HermMatrix<T,U,S,I>& operator=(T x) 
    {
      TMVAssert(IMAG(x) == RealType(T)(0));
      return SetToIdentity(x); 
    }

    inline HermMatrix<T,U,S,I>& operator=(
        const AssignableToSymMatrix<RealType(T)>& m2) 
    { 
      TMVAssert(size() == m2.colsize());
      TMVAssert(m2.isherm());
      m2.AssignToS(View());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    inline HermMatrix<T,U,S,I>& operator=(
        const AssignableToSymMatrix<ComplexType(T)>& m2) 
    { 
      TMVAssert(IsComplex(T()));
      TMVAssert(size() == m2.colsize());
      TMVAssert(m2.isherm());
      m2.AssignToS(View());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    inline HermMatrix<T,U,S,I>& operator=(
        const GenDiagMatrix<RealType(T)>& m2) 
    { 
      TMVAssert(size() == m2.colsize());
      View() = m2;
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    inline HermMatrix<T,U,S,I>& operator=(
        const GenDiagMatrix<ComplexType(T)>& m2) 
    { TMVAssert(FALSE); return *this; }

    template <class T2> 
    inline HermMatrix<T,U,S,I>& operator=(const OProdVV<T,T2,T2>& opvv)
    {
      TMVAssert(size() == opvv.colsize());
      TMVAssert(size() == opvv.rowsize());
      TMVAssert(opvv.GetV1().SameAs(opvv.GetV2().Conjugate()));
      Rank1Update<false>(opvv.GetX(),opvv.GetV1(),View());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    template <class T2> 
    inline HermMatrix<T,U,S,I>& operator=(const ProdMM<T,T2,T2>& pmm) 
    {
      TMVAssert(size() == pmm.colsize());
      TMVAssert(size() == pmm.rowsize());
      TMVAssert(pmm.GetM1().SameAs(pmm.GetM2().Adjoint()));
      RankKUpdate<false>(pmm.GetX(),pmm.GetM1(),View());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    template <class T2> 
    inline HermMatrix<T,U,S,I>& operator=(const ProdUL<T,T2,T2>& pmm)
    {
      TMVAssert(size() == pmm.colsize());
      TMVAssert(size() == pmm.rowsize());
      TMVAssert(pmm.GetM1().SameAs(pmm.GetM2().Adjoint()));
      RankKUpdate<false>(pmm.GetX(),pmm.GetM1(),View());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    template <class T2> 
    inline HermMatrix<T,U,S,I>& operator=(const ProdLU<T,T2,T2>& pmm) 
    {
      TMVAssert(size() == pmm.colsize());
      TMVAssert(size() == pmm.rowsize());
      TMVAssert(pmm.GetM1().SameAs(pmm.GetM2().Adjoint()));
      RankKUpdate<false>(pmm.GetX(),pmm.GetM1(),View());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }


    //
    // Access
    //

    typedef RefType(T) reference;

    inline T operator()(int i, int j) const
    {
      if (I==CStyle) {
        TMVAssert(i>=0 && i<int(size()));
        TMVAssert(j>=0 && j<int(size()));
        return cref(i,j); 
      } else {
        TMVAssert(i>0 && i<=int(size()));
        TMVAssert(j>0 && j<=int(size()));
        return cref(i-1,j-1); 
      }
    }

    inline reference operator()(int i, int j) 
    { 
      if (I==CStyle) {
        TMVAssert(i>=0 && i<int(size()));
        TMVAssert(j>=0 && j<int(size()));
        return ref(i,j); 
      } else {
        TMVAssert(i>0 && i<=int(size()));
        TMVAssert(j>0 && j<=int(size()));
        return ref(i-1,j-1); 
      }
    }

    inline ConstVectorView<T,I> row(int i, int j1, int j2) const 
    { 
      if (I==FortranStyle) {
        TMVAssert(i>0 && i<=int(size())); --i;
        TMVAssert(j1>0 && j1-j2<=0 && j2<=int(size())); --j1;
      } else {
        TMVAssert(i>=0 && i<int(size()));
        TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size()));
      }
      TMVAssert( i-j1<=0 || j2-i<=1 );
      if ((U==Upper && i-j1<=0) || (U==Lower && j2-i<=1))
        return ConstVectorView<T,I>(itsm.get()+i*stepi()+j1*stepj(),
            j2-j1,stepj(),NonConj); 
      else
        return ConstVectorView<T,I>(itsm.get()+i*stepj()+j1*stepi(),
            j2-j1,stepi(),ConjOf(T,NonConj)); 
    }

    inline ConstVectorView<T,I> col(int j, int i1, int i2) const
    {
      if (I==FortranStyle) { 
        TMVAssert(j>0 && j<=int(size())); --j;
        TMVAssert(i1>0 && i1-i2<=0 && i2<=int(size())); --i1;
      } else {
        TMVAssert(j>=0 && j<int(size()));
        TMVAssert(i1>=0 && i1-i2<=0 && i2<=int(size()));
      }
      TMVAssert( j-i1<=0 || i2-j<=1 );
      if ((U==Upper && i2-j<=1) || (U==Lower && j-i1<=0))
        return ConstVectorView<T,I>(itsm.get()+i1*stepi()+j*stepj(),
            i2-i1,stepi(),NonConj); 
      else 
        return ConstVectorView<T,I>(itsm.get()+i1*stepj()+j*stepi(),
            i2-i1,stepj(),ConjOf(T,NonConj)); 
    }

    inline ConstVectorView<T,I> diag() const
    { return ConstVectorView<T,I>(itsm.get(),size(),stepi()+stepj(),NonConj); }

    inline ConstVectorView<T,I> diag(int i) const
    {
      TMVAssert(i>=-int(size()) && i<=int(size())); 
      ConjType newct = 
      ((i>0) == (U==Upper)) ? NonConj : ConjOf(T,NonConj);
      i = std::abs(i);
      TMVAssert(i<=int(size())); 
      return ConstVectorView<T,I>(itsm.get()+i*(U==Upper?stepj():stepi()),
          size()-i,stepi()+stepj(),newct);
    }

    inline ConstVectorView<T,I> diag(int i, int j1, int j2) const
    {
      TMVAssert(i>=-int(size()) && i<=int(size())); 
      if (I==FortranStyle) {
        TMVAssert(j1>0 && j1-j2<=0 && j2<=int(size())-std::abs(i)); --j1;
      } else {
        TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size())-std::abs(i));
      }
      ConjType newct = 
      ((i>0) == (U==Upper)) ? NonConj : ConjOf(T,NonConj);
      i = std::abs(i);
      const int ds = stepi()+stepj();
      return ConstVectorView<T,I>(itsm.get()+i*(U==Upper?stepj():stepi())+j1*ds,
          j2-j1,ds,newct);
    }

    inline VectorView<T,I> row(int i, int j1, int j2)
    { 
      if (I==FortranStyle) {
        TMVAssert(i>0 && i<=int(size())); --i;
        TMVAssert(j1>0 && j1-j2<=0 && j2<=int(size())); --j1;
      } else {
        TMVAssert(i>=0 && i<int(size()));
        TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size()));
      }
      TMVAssert( i-j1<=0 || j2-i<=1 );
      if ((U==Upper && i-j1<=0) || (U==Lower && j2-i<=1))
        return VectorView<T,I>(itsm.get()+i*stepi()+j1*stepj(),
            j2-j1,stepj(),NonConj FIRSTLAST); 
      else
        return VectorView<T,I>(itsm.get()+i*stepj()+j1*stepi(),
            j2-j1,stepi(),ConjOf(T,NonConj) FIRSTLAST); 
    }

    inline VectorView<T,I> col(int j, int i1, int i2)
    {
      if (I==FortranStyle) { 
        TMVAssert(j>0 && j<=int(size())); --j;
        TMVAssert(i1>0 && i1-i2<=0 && i2<=int(size())); --i1;
      } else {
        TMVAssert(j>=0 && j<int(size()));
        TMVAssert(i1>=0 && i1-i2<=0 && i2<=int(size()));
      }
      TMVAssert( j-i1<=0 || i2-j<=1 );
      if ((U==Upper && i2-j<=1) || (U==Lower && j-i1<=0))
        return VectorView<T,I>(itsm.get()+i1*stepi()+j*stepj(),
            i2-i1,stepi(),NonConj FIRSTLAST); 
      else 
        return VectorView<T,I>(itsm.get()+i1*stepj()+j*stepi(),
            i2-i1,stepj(),ConjOf(T,NonConj) FIRSTLAST); 
    }

    inline VectorView<T,I> diag()
    {
      return VectorView<T,I>(itsm.get(),size(),stepi()+stepj(),NonConj FIRSTLAST);
    }

    inline VectorView<T,I> diag(int i) 
    {
      TMVAssert(i>=-int(size()) && i<=int(size())); 
      ConjType newct = 
      ((i>0) == (U==Upper)) ? NonConj : ConjOf(T,NonConj);
      i = std::abs(i);
      TMVAssert(i<=int(size())); 
      return VectorView<T,I>(itsm.get()+i*(U==Upper?stepj():stepi()),size()-i,
          stepi()+stepj(),newct FIRSTLAST); 
    }

    inline VectorView<T,I> diag(int i, int j1, int j2) 
    {
      TMVAssert(i>=-int(size()) && i<=int(size())); 
      if (I==FortranStyle) {
        TMVAssert(j1>0 && j1-j2<=0 && j2<=int(size())-std::abs(i)); --j1;
      } else {
        TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size())-std::abs(i));
      }
      ConjType newct = 
      ((i>0) == (U==Upper)) ? NonConj : ConjOf(T,NonConj);
      i = std::abs(i);
      const int ds = stepi()+stepj();
      return VectorView<T,I>(itsm.get()+i*(U==Upper?stepj():stepi())+j1*ds,
          j2-j1,ds,newct FIRSTLAST);
    }

    //
    // Modifying Functions
    //

    inline HermMatrix<T,U,S,I>& Zero() 
    {
      std::fill_n(itsm.get(),itslen,T(0));
      return *this;
    }

    inline HermMatrix<T,U,S,I>& SetAllTo(T x) 
    { 
      TMVAssert(IMAG(x) == RealType(T)(0));
      UpperTri().SetAllTo(x); 
      return *this; 
    }

    inline HermMatrix<T,U,S,I>& Clip(RealType(T) thresh) 
    { UpperTri().Clip(thresh); return *this; }

    inline HermMatrix<T,U,S,I>& ConjugateSelf() 
    { if (IsComplex(T())) UpperTri().ConjugateSelf(); return *this; }

    inline HermMatrix<T,U,S,I>& TransposeSelf() 
    { if (IsComplex(T())) UpperTri().ConjugateSelf(); return *this; }

    inline HermMatrix<T,U,S,I>& SetToIdentity(T x=T(1)) 
    { 
      TMVAssert(IMAG(x) == RealType(T)(0));
      Zero(); 
      diag().SetAllTo(x); 
      return *this; 
    }

    inline HermMatrix<T,U,S,I>& SwapRowsCols(int i1, int i2)
    { View().SwapRowsCols(i1,i2); return *this; }

    inline HermMatrix<T,U,S,I>& PermuteRowsCols(
        const int* p, int i1, int i2)
    { View().PermuteRowsCols(p,i1,i2); return *this; }

    inline HermMatrix<T,U,S,I>& ReversePermuteRowsCols(
        const int* p, int i1, int i2)
    { View().ReversePermuteRowsCols(p,i1,i2); return *this; }

    inline HermMatrix<T,U,S,I>& PermuteRowsCols(const int* p)
    { View().PermuteRowsCols(p); return *this; }

    inline HermMatrix<T,U,S,I>& ReversePermuteRowsCols(const int* p)
    { View().ReversePermuteRowsCols(p); return *this; }


    //
    // SubMatrix
    //

    inline ConstMatrixView<T,I> SubMatrix(
        int i1, int i2, int j1, int j2) const
    {
      TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,1,1));
      if (I==FortranStyle) { --i1; --j1; }
      TMVAssert(i2-1<=j1 || j2-1<=i1);
      if ((U==Upper && i2-1<=j1) || (U==Lower && j2-1<=i1))
        return ConstMatrixView<T,I>(itsm.get()+i1*stepi()+j1*stepj(),
            i2-i1,j2-j1,stepi(),stepj(),S,NonConj);
      else
        return ConstMatrixView<T,I>(itsm.get()+i1*stepj()+j1*stepi(),
            i2-i1,j2-j1,stepj(),stepi(),TransOf(S),ConjOf(T,NonConj));
    }

    inline ConstMatrixView<T,I> SubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) const
    {
      const StorageType newstor = S==RowMajor ?
      jstep == 1 ? RowMajor : NoMajor :
      istep == 1 ? ColMajor : NoMajor;
      TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,istep,jstep));
      if (I==FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
      TMVAssert(i2-istep<=j1 || j2-jstep<=i1);
      if ((U==Upper && i2-istep<=j1) || (U==Lower && j2-jstep<=i1))
        return ConstMatrixView<T,I>(itsm.get()+i1*stepi()+j1*stepj(),
            (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
            newstor,NonConj);
      else
        return ConstMatrixView<T,I>(itsm.get()+i1*stepj()+j1*stepi(),
            (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
            TransOf(newstor),ConjOf(T,NonConj));
    }

    inline ConstVectorView<T,I> SubVector(
        int i, int j, int istep, int jstep, int n) const
    {
      TMVAssert(GenSymMatrix<T>::OKSubVector(i,j,istep,jstep,n));
      if ((U==Upper && i-j<=0) || (U==Lower && j-i<=0))
        return ConstVectorView<T,I>(itsm.get()+i*stepi()+j*stepj(),n,
            istep*stepi()+jstep*stepj(),NonConj);
      else
        return ConstVectorView<T,I>(itsm.get()+i*stepj()+j*stepi(),n,
            istep*stepj()+jstep*stepi(),ConjOf(T,NonConj));
    }

    inline ConstSymMatrixView<T,I> SubSymMatrix(int i1, int i2) const
    {
      TMVAssert(View().OKSubSymMatrix(i1,i2,1));
      if (I==FortranStyle) { --i1; }
      return ConstSymMatrixView<T,I>(itsm.get()+i1*(stepi()+stepj()),
          i2-i1,stepi(),stepj(),Herm,U,S,NonConj);
    }

    inline ConstSymMatrixView<T,I> SubSymMatrix(
        int i1, int i2, int istep) const
    {
      TMVAssert(View().OKSubSymMatrix(i1,i2,istep));
      if (I==FortranStyle) { --i1; i2+=istep-1; }
      return ConstSymMatrixView<T,I>(itsm.get()+i1*(stepi()+stepj()),
          (i2-i1)/istep,istep*stepi(),istep*stepj(),Herm,U,
          istep==1 ? S : NoMajor, NonConj);
    }

    inline ConstUpperTriMatrixView<T,I> UpperTri(
        DiagType dt = NonUnitDiag) const
    {
      return U==Upper ? 
      ConstUpperTriMatrixView<T,I>(itsm.get(),size(),stepi(),stepj(),
          dt,S,NonConj) :
      ConstUpperTriMatrixView<T,I>(itsm.get(),size(),stepj(),stepi(),
          dt,TransOf(S),ConjOf(T,NonConj));
    }

    inline ConstLowerTriMatrixView<T,I> LowerTri(
        DiagType dt = NonUnitDiag) const
    {
      return U==Lower ? 
      ConstLowerTriMatrixView<T,I>(itsm.get(),size(),stepi(),stepj(),
          dt,S,NonConj) :
      ConstLowerTriMatrixView<T,I>(itsm.get(),size(),stepj(),stepi(),
          dt,TransOf(S),ConjOf(T,NonConj));
    }

    inline MatrixView<T,I> SubMatrix(int i1, int i2, int j1, int j2)
    {
      TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,1,1));
      if (I==FortranStyle) { --i1; --j1; }
      TMVAssert(i2-1<=j1 || j2-1<=i1);
      if ((U==Upper && i2-1<=j1) || (U==Lower && j2-1<=i1))
        return MatrixView<T,I>(itsm.get()+i1*stepi()+j1*stepj(),
            i2-i1,j2-j1,stepi(),stepj(),S,NonConj FIRSTLAST);
      else
        return MatrixView<T,I>(itsm.get()+i1*stepj()+j1*stepi(),
            i2-i1,j2-j1,stepj(),stepi(),TransOf(S),ConjOf(T,NonConj) 
            FIRSTLAST);
    }

    inline MatrixView<T,I> SubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep)
    {
      const StorageType newstor = S==RowMajor ?
      jstep == 1 ? RowMajor : NoMajor :
      istep == 1 ? ColMajor : NoMajor;
      TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,istep,jstep));
      if (I==FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
      TMVAssert(i2-istep<=j1 || j2-jstep<=i1);
      if ((U==Upper && i2-istep<=j1) || (U==Lower && j2-jstep<=i1))
        return MatrixView<T,I>(itsm.get()+i1*stepi()+j1*stepj(),
            (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
            newstor,NonConj FIRSTLAST);
      else
        return MatrixView<T,I>(itsm.get()+i1*stepj()+j1*stepi(),
            (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
            TransOf(newstor),ConjOf(T,NonConj) FIRSTLAST);
    }

    inline VectorView<T,I> SubVector(
        int i, int j, int istep, int jstep, int n)
    {
      TMVAssert(GenSymMatrix<T>::OKSubVector(i,j,istep,jstep,n));
      if ((U==Upper && i-j<=0) || (U==Lower && j-i<=0))
        return VectorView<T,I>(itsm.get()+i*stepi()+j*stepj(),n,
            istep*stepi()+jstep*stepj(),NonConj FIRSTLAST);
      else
        return VectorView<T,I>(itsm.get()+i*stepj()+j*stepi(),n,
            istep*stepj()+jstep*stepi(),ConjOf(T,NonConj) FIRSTLAST);
    }

    inline SymMatrixView<T,I> SubSymMatrix(int i1, int i2)
    {
      TMVAssert(View().OKSubSymMatrix(i1,i2,1));
      if (I==FortranStyle) { --i1; }
      return SymMatrixView<T,I>(itsm.get()+i1*(stepi()+stepj()),
          i2-i1,stepi(),stepj(),Herm,U,S,NonConj FIRSTLAST);
    }

    inline SymMatrixView<T,I> SubSymMatrix(int i1, int i2, int istep) 
    {
      TMVAssert(View().OKSubSymMatrix(i1,i2,istep));
      if (I==FortranStyle) { --i1; i2+=istep-1; }
      return SymMatrixView<T,I>(itsm.get()+i1*(stepi()+stepj()),
          (i2-i1)/istep,istep*stepi(),istep*stepj(),Herm,U,
          istep==1 ? S : NoMajor,NonConj FIRSTLAST);
    }

    inline UpperTriMatrixView<T,I> UpperTri(DiagType dt = NonUnitDiag)
    {
      return U==Upper ? 
      UpperTriMatrixView<T,I>(itsm.get(),size(),stepi(),stepj(),
          dt,S,NonConj FIRSTLAST) :
      UpperTriMatrixView<T,I>(itsm.get(),size(),stepj(),stepi(),
          dt,TransOf(S),ConjOf(T,NonConj) FIRSTLAST);
    }

    inline LowerTriMatrixView<T,I> LowerTri(DiagType dt = NonUnitDiag)
    {
      return U==Lower ? 
      LowerTriMatrixView<T,I>(itsm.get(),size(),stepi(),stepj(),
          dt,S,NonConj FIRSTLAST) :
      LowerTriMatrixView<T,I>(itsm.get(),size(),stepj(),stepi(),
          dt,TransOf(S),ConjOf(T,NonConj) FIRSTLAST);
    }

    inline ConstSymMatrixView<T,I> View() const
    { 
      return ConstSymMatrixView<T,I>(itsm.get(),size(),stepi(),stepj(),
          Herm,U,S,NonConj);
    }

    inline ConstSymMatrixView<T,I> Transpose() const
    {
      return ConstSymMatrixView<T,I>(itsm.get(),size(),stepj(),stepi(),
          Herm,UTransOf(U),TransOf(S),NonConj);
    }

    inline ConstSymMatrixView<T,I> Conjugate() const
    { 
      return ConstSymMatrixView<T,I>(itsm.get(),size(),stepi(),stepj(),
          Herm,U,S,ConjOf(T,NonConj));
    }

    inline ConstSymMatrixView<T,I> Adjoint() const
    {
      return ConstSymMatrixView<T,I>(itsm.get(),size(),stepj(),stepi(),
          Herm,UTransOf(U),TransOf(S),ConjOf(T,NonConj));
    }

    inline ConstSymMatrixView<RealType(T),I> Real() const
    {
      return ConstSymMatrixView<RealType(T),I>(
          reinterpret_cast<const RealType(T)*>(itsm.get()),size(),
          IsReal(T()) ? stepi() : 2*stepi(),
          IsReal(T()) ? stepj() : 2*stepj(),
          Herm,U,IsReal(T())?S:NoMajor,NonConj);
    }

    inline ConstSymMatrixView<RealType(T),I> Imag() const
    // The imaginary part of a Hermitian matrix is anti-symmetric
    // so this is illegal.
    { 
      TMVAssert(FALSE);
      return ConstSymMatrixView<RealType(T),I>(0,0,0,0,Herm,U,S,NonConj);
    }

    inline SymMatrixView<T,I> View() 
    { 
      return SymMatrixView<T,I>(itsm.get(),size(),stepi(),stepj(),
          Herm,U,S,NonConj FIRSTLAST);
    }

    inline SymMatrixView<T,I> Transpose() 
    {
      return SymMatrixView<T,I>(itsm.get(),size(),stepj(),stepi(),
          Herm,UTransOf(U),TransOf(S),NonConj FIRSTLAST);
    }

    inline SymMatrixView<T,I> Conjugate() 
    { 
      return SymMatrixView<T,I>(itsm.get(),size(),stepi(),stepj(),
          Herm,U,S,ConjOf(T,NonConj) FIRSTLAST);
    }

    inline SymMatrixView<T,I> Adjoint() 
    {
      return SymMatrixView<T,I>(itsm.get(),size(),stepj(),stepi(),
          Herm,UTransOf(U),TransOf(S),ConjOf(T,NonConj) FIRSTLAST);
    }

    inline SymMatrixView<RealType(T),I> Real()
    {
      return SymMatrixView<RealType(T),I>(
          reinterpret_cast<RealType(T)*>(itsm.get()),size(),
          IsReal(T()) ? stepi() : 2*stepi(),
          IsReal(T()) ? stepj() : 2*stepj(),
          Herm,U,IsReal(T())?S:NoMajor,NonConj
#ifdef TMVFLDEBUG
          ,reinterpret_cast<const RealType(T)*>(first)
          ,reinterpret_cast<const RealType(T)*>(last)
#endif
          );
    }

    inline SymMatrixView<RealType(T),I> Imag()
    // The imaginary part of a Hermitian matrix is anti-symmetric
    // so this is illegal.
    { 
      TMVAssert(FALSE);
      return SymMatrixView<RealType(T),I>(0,0,0,0,Herm,U,S,NonConj
          FIRSTLAST1(0,0) );
    }

    inline size_t size() const { return itss; }
    inline const T* cptr() const { return itsm.get(); }
    inline T* ptr() { return itsm.get(); }
    inline int stepi() const { return S==RowMajor ? itss : 1; }
    inline int stepj() const { return S==RowMajor ? 1 : itss; }
    inline SymType sym() const { return Herm; }
    inline UpLoType uplo() const { return U; }
    inline StorageType stor() const { return S; }
    inline ConjType ct() const { return NonConj; }
    inline bool isrm() const { return S==RowMajor; }
    inline bool iscm() const { return S==ColMajor; }
    inline bool isconj() const { return false; }
    inline bool isherm() const { return true; }
    inline bool issym() const { return IsReal(T()); }
    inline bool isupper() const { return U == Upper; }

    inline reference ref(int i, int j)
    {
      if ((U==Upper && i <= j) || (U==Lower && i>=j)) 
        if (S == RowMajor) return REF(itsm.get() + i*itss + j, NonConj);
        else return REF(itsm.get() + j*itss + i, NonConj);
      else 
        if (S == RowMajor) return REF(itsm.get() + j*itss + i, Conj);
        else return REF(itsm.get() + i*itss + j, Conj);
    }

    inline T cref(int i, int j) const 
    {
      if ((U==Upper && i <= j) || (U==Lower && i>=j)) 
        if (S == RowMajor) return itsm.get()[i*itss + j];
        else return itsm.get()[j*itss + i];
      else 
        if (S == RowMajor) return CONJ(itsm.get()[j*itss + i]);
        else return CONJ(itsm.get()[i*itss + j]);
    }

  protected :

    const size_t itslen;
    auto_array<T> itsm;
    const size_t itss;

#ifdef TMVFLDEBUG
  public:
    const T*const first;
    const T*const last;
#endif

  }; // HermMatrix

  //---------------------------------------------------------------------------

  //
  // Special Creators: 
  //   SymMatrixViewOf(m,uplo)
  //   HermMatrixViewOf(m,uplo)
  //   SymMatrixViewOf(mptr,size,uplo,stor)
  //   HermMatrixViewOf(mptr,size,uplo,stor)
  //

  template <class T> 
  inline ConstSymMatrixView<T> SymMatrixViewOf(
      const GenMatrix<T>& m, UpLoType uplo)
  {
    TMVAssert(m.colsize()==m.rowsize());
    return ConstSymMatrixView<T>(m.cptr(),m.rowsize(),m.stepi(),m.stepj(),
        Sym,uplo,m.stor(),m.ct()); 
  }

  template <class T, IndexStyle I> 
  inline ConstSymMatrixView<T,I> SymMatrixViewOf(
      const ConstMatrixView<T,I>& m, UpLoType uplo)
  { 
    TMVAssert(m.colsize()==m.rowsize());
    return ConstSymMatrixView<T,I>(m.cptr(),m.rowsize(),m.stepi(),m.stepj(),
        Sym,uplo,m.stor(),m.ct()); 
  }

  template <class T, StorageType S, IndexStyle I> 
  inline ConstSymMatrixView<T,I> SymMatrixViewOf(
      const Matrix<T,S,I>& m, UpLoType uplo)
  {
    TMVAssert(m.colsize()==m.rowsize());
    return ConstSymMatrixView<T,I>(m.cptr(),m.rowsize(),m.stepi(),m.stepj(),
        Sym,uplo,m.stor(),m.ct()); 
  }

  template <class T, IndexStyle I> 
  inline SymMatrixView<T,I> SymMatrixViewOf(
      const MatrixView<T,I>& m, UpLoType uplo)
  { 
    TMVAssert(m.colsize()==m.rowsize());
    return SymMatrixView<T,I>(m.ptr(),m.rowsize(),m.stepi(),m.stepj(),
        Sym,uplo,m.stor(),m.ct() FIRSTLAST1(m.first,m.last)); 
  }

  template <class T, StorageType S, IndexStyle I> 
  inline SymMatrixView<T,I> SymMatrixViewOf(
      Matrix<T,S,I>& m, UpLoType uplo)
  {
    TMVAssert(m.colsize()==m.rowsize());
    return SymMatrixView<T,I>(m.ptr(),m.rowsize(),m.stepi(),m.stepj(),
        Sym,uplo,m.stor(),m.ct() FIRSTLAST1(m.first,m.last)); 
  }

  template <class T> 
  inline ConstSymMatrixView<T> HermMatrixViewOf(
      const GenMatrix<T>& m, UpLoType uplo)
  {
    TMVAssert(m.colsize()==m.rowsize());
#ifdef XTEST
    TMVAssert(IsReal(T()) || m.diag().Imag().NormInf() == RealType(T)(0));
#endif
    return ConstSymMatrixView<T>(m.cptr(),m.rowsize(),m.stepi(),m.stepj(),
        Herm,uplo,m.stor(),m.ct()); 
  }

  template <class T, IndexStyle I> 
  inline ConstSymMatrixView<T> HermMatrixViewOf(
      const ConstMatrixView<T,I>& m, UpLoType uplo)
  { 
    TMVAssert(m.colsize()==m.rowsize());
#ifdef XTEST
    TMVAssert(IsReal(T()) || m.diag().Imag().NormInf() == RealType(T)(0));
#endif
    return ConstSymMatrixView<T,I>(m.cptr(),m.rowsize(),m.stepi(),m.stepj(),
        Herm,uplo,m.stor(),m.ct()); 
  }

  template <class T, StorageType S, IndexStyle I> 
  inline ConstSymMatrixView<T,I> HermMatrixViewOf(
      const Matrix<T,S,I>& m, UpLoType uplo)
  {
    TMVAssert(m.colsize()==m.rowsize());
#ifdef XTEST
    TMVAssert(IsReal(T()) || m.diag().Imag().NormInf() == RealType(T)(0));
#endif
    return ConstSymMatrixView<T,I>(m.cptr(),m.rowsize(),m.stepi(),m.stepj(),
        Herm,uplo,m.stor(),m.ct()); 
  }

  template <class T, IndexStyle I> 
  inline SymMatrixView<T> HermMatrixViewOf(
      const MatrixView<T,I>& m, UpLoType uplo)
  { 
    TMVAssert(m.colsize()==m.rowsize());
#ifdef XTEST
    TMVAssert(IsReal(T()) || m.diag().Imag().NormInf() == RealType(T)(0));
#endif
    return SymMatrixView<T,I>(m.ptr(),m.rowsize(),m.stepi(),m.stepj(),
        Herm,uplo,m.stor(),m.ct() FIRSTLAST1(m.first,m.last)); 
  }

  template <class T, StorageType S, IndexStyle I> 
  inline SymMatrixView<T,I> HermMatrixViewOf(
      Matrix<T,S,I>& m, UpLoType uplo)
  {
    TMVAssert(m.colsize()==m.rowsize());
#ifdef XTEST
    TMVAssert(IsReal(T()) || m.diag().Imag().NormInf() == RealType(T)(0));
#endif
    return SymMatrixView<T,I>(m.ptr(),m.rowsize(),m.stepi(),m.stepj(),
        Herm,uplo,m.stor(),m.ct() FIRSTLAST1(m.first,m.last)); 
  }

  template <class T> 
  ConstSymMatrixView<T> SymMatrixViewOf(
      const T* vv, size_t size, UpLoType uplo, StorageType stor);

  template <class T> 
  ConstSymMatrixView<T> HermMatrixViewOf(
      const T* vv, size_t size, UpLoType uplo, StorageType stor);

  template <class T> 
  SymMatrixView<T> SymMatrixViewOf(
      T* vv, size_t size, UpLoType uplo, StorageType stor);

  template <class T> 
  SymMatrixView<T> HermMatrixViewOf(
      T* vv, size_t size, UpLoType uplo, StorageType stor);

  //
  // Swap Matrices
  //

  template <class T> 
  inline void Swap(const SymMatrixView<T>& m1, const SymMatrixView<T>& m2)
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(m1.issym() == m2.issym());
    TMVAssert(m1.isherm() == m2.isherm());
    Swap(m1.UpperTri(),m2.UpperTri()); 
  }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline void Swap(const SymMatrixView<T>& m1, SymMatrix<T,U,S,I>& m2)
  { Swap(m1,m2.View()); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline void Swap(SymMatrix<T,U,S,I>& m1, const SymMatrixView<T>& m2)
  { Swap(m1.View(),m2); }

  template <class T, UpLoType U1, StorageType S1, IndexStyle I1, UpLoType U2, StorageType S2, IndexStyle I2>
  inline void Swap(SymMatrix<T,U1,S1,I1>& m1, SymMatrix<T,U2,S2,I2>& m2)
  { Swap(m1.View(),m2.View()); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline void Swap(const SymMatrixView<T>& m1, HermMatrix<T,U,S,I>& m2)
  { Swap(m1,m2.View()); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline void Swap(HermMatrix<T,U,S,I>& m1, const SymMatrixView<T>& m2)
  { Swap(m1.View(),m2); }

  template <class T, UpLoType U1, StorageType S1, IndexStyle I1, UpLoType U2, StorageType S2, IndexStyle I2>
  inline void Swap(HermMatrix<T,U1,S1,I1>& m1, HermMatrix<T,U2,S2,I2>& m2)
  { Swap(m1.View(),m2.View()); }


  //
  // Views:
  //

  template <class T> 
  inline ConstSymMatrixView<T> Transpose(const GenSymMatrix<T>& m)
  { return m.Transpose(); }

  template <class T, IndexStyle I> 
  inline ConstSymMatrixView<T,I> Transpose(const ConstSymMatrixView<T,I>& m)
  { return m.Transpose(); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline ConstSymMatrixView<T,I> Transpose(const SymMatrix<T,U,S,I>& m)
  { return m.Transpose(); }

  template <class T, IndexStyle I> 
  inline SymMatrixView<T,I> Transpose(const SymMatrixView<T,I>& m)
  { return m.Transpose(); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline SymMatrixView<T,I> Transpose(SymMatrix<T,U,S,I>& m)
  { return m.Transpose(); }

  template <class T> 
  inline ConstSymMatrixView<T> Conjugate(const GenSymMatrix<T>& m)
  { return m.Conjugate(); }

  template <class T, IndexStyle I> 
  inline ConstSymMatrixView<T,I> Conjugate(const ConstSymMatrixView<T,I>& m)
  { return m.Conjugate(); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline ConstSymMatrixView<T,I> Conjugate(const SymMatrix<T,U,S,I>& m)
  { return m.Conjugate(); }

  template <class T, IndexStyle I> 
  inline SymMatrixView<T,I> Conjugate(const SymMatrixView<T,I>& m)
  { return m.Conjugate(); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline SymMatrixView<T,I> Conjugate(SymMatrix<T,U,S,I>& m)
  { return m.Conjugate(); }

  template <class T> 
  inline ConstSymMatrixView<T> Adjoint(const GenSymMatrix<T>& m)
  { return m.Adjoint(); }

  template <class T, IndexStyle I> 
  inline ConstSymMatrixView<T,I> Adjoint(const ConstSymMatrixView<T,I>& m)
  { return m.Adjoint(); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline ConstSymMatrixView<T,I> Adjoint(const SymMatrix<T,U,S,I>& m)
  { return m.Adjoint(); }

  template <class T, IndexStyle I> 
  inline SymMatrixView<T,I> Adjoint(const SymMatrixView<T,I>& m)
  { return m.Adjoint(); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline SymMatrixView<T,I> Adjoint(SymMatrix<T,U,S,I>& m)
  { return m.Adjoint(); }

  template <class T> 
  inline QuotXS<T,T> Inverse(const GenSymMatrix<T>& m)
  { return m.Inverse(); }

  //
  // SymMatrix ==, != SymMatrix
  //

  template <class T1, class T2> 
  inline bool operator==(
      const GenSymMatrix<T1>& m1, const GenSymMatrix<T2>& m2)
  { return m1.UpperTri() == m2.UpperTri(); }
  template <class T1, class T2> 
  inline bool operator!=(
      const GenSymMatrix<T1>& m1, const GenSymMatrix<T2>& m2)
  { return !(m1 == m2); }


  //
  // I/O
  //

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  std::istream& operator>>(std::istream& is, 
      auto_ptr<SymMatrix<T,U,S,I> >& m);

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  std::istream& operator>>(std::istream& is, 
      auto_ptr<HermMatrix<T,U,S,I> >& m);

  template <class T> 
  std::istream& operator>>(std::istream& is, const SymMatrixView<T>& m);

  template <class T, UpLoType U, StorageType S, IndexStyle I>
  inline std::istream& operator>>(std::istream& is, SymMatrix<T,U,S,I>& m)
  { return is>>m.View(); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline std::istream& operator>>(std::istream& is, HermMatrix<T,U,S,I>& m)
  { return is>>m.View(); }

} // namespace tmv

#endif
