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
// This file defines the TMV SymBandMatrix and HermBandMatrix classes.
//
// Constructors:
//
//    SymBandMatrix is used for symmetric band matrices, and 
//    HermBandMatrix is used for Hermitian band matrices.  
//    For real matrices, these are the same thing:
//        A = A.Transpose().
//    But for complex, they are different:
//        A_sym = A_sym.Transpose()
//        A_herm = A_herm.Adjoint()
//
//    For these notes, I will always write SymBandMatrix, but (except where
//    otherwise indicated) everything applies the same for Sym and Herm.
//    Also, the Views keep track of sym/herm difference with a parameter,
//    so it is always a GenSymBandMatrix, ConstSymBandMatrixView, or 
//    SymBandMatrixView - never Herm in any of these.
//
//    Caveat: Complex Hermitian matrices are such that A = At, which 
//    implies that their diagonal elements are real.  Many routines
//    involving HermBandMatrixes assume the reality of the diagonal.
//    However, it is possible to assign a non-real value to a diagonal
//    element.  If the user does this, the results are undefined.
//
//    In addition to the type template parameter (T), SymBandMatrixes have two
//    additional template parameters:
//        UpLoType uplo = Upper || Lower
//        StorageType stor = RowMajor || ColMajor || DiagMajor
//
//        They both have default values, so you can omit both, or
//        just stor.  The default values are: {Upper, RowMajor}
//
//        The first, uplo, refers to which triangular half stores the actual 
//        data.
//
//        The storage follows the same meaning as for regular Matrices.
//
//    SymBandMatrix<T,uplo,stor>(size_t n, int nlo)
//        Makes a SymBandMatrix with column size = row size = n
//        and nlo (=nhi) off-diagonals with _uninitialized_ values.
//
//    SymBandMatrix<T,uplo,stor>(size_t n, int nlo, T x)
//        Makes a SymBandMatrix with values all initialized to x
//        For Hermitian matrixces, x must be real.
//
//    SymBandMatrix<T,uplo,stor>(size_t n, int nlo, const T* m)
//    SymBandMatrix<T,uplo,stor>(size_t n, int nlo, const vector<T>& m)
//        Make a SymBandMatrix with copies the elements of m which
//        fall in tha appropriate upper or lower triangle.
//        The lengths of the arrays must be BandStorageLength(stor,n,n,nlo,0).
//
//    SymBandMatrix<T,uplo,stor>(const Matrix<T>& m)
//    SymBandMatrix<T,uplo,stor>(const SymMatrix<T>& m)
//    SymBandMatrix<T,uplo,stor>(const BandMatrix<T>& m)
//        Makes a SymBandMatrix which copies the corresponding elements of m.
//
//    ConstSymBandMatrixView<T> SymBandMatrixViewOf(const Matrix<T>& m, 
//        uplo, int nlo)
//    ConstSymBandMatrixView<T> HermBandMatrixViewOf(const Matrix<T>& m, 
//        uplo, int nlo)
//    ConstSymBandMatrixView<T> SymBandMatrixViewOf(const BandMatrix<T>& m, 
//        uplo)
//    ConstSymBandMatrixView<T> HermBandMatrixViewOf(const BandMatrix<T>& m, 
//        uplo)
//    ConstSymBandMatrixView<T> SymBandMatrixViewOf(const BandMatrix<T>& m, 
//        uplo, int nlo)
//    ConstSymBandMatrixView<T> HermBandMatrixViewOf(const BandMatrix<T>& m, 
//        uplo, int nlo)
//    ConstSymBandMatrixView<T> SymBandMatrixViewOf(const SymMatrix<T>& m, 
//        int nlo)
//    ConstSymBandMatrixView<T> SymBandMatrixViewOf(const HermMatrix<T>& m, 
//        int nlo)
//        Makes a constant SymBandMatrix view of the corresponding part of m.
//        While this view cannot be modified, changing the original m
//        will cause corresponding changes in this view of m.
//
//    SymMatrixView<T> SymMatrixViewOf(Matrix<T>& m, uplo, int nlo)
//    SymMatrixView<T> HermMatrixViewOf(Matrix<T>& m, uplo, int nlo)
//    [ ... ]
//        Makes a modifiable SymMatrix view of the corresponding part of m.
//        Modifying this matrix will change the corresponding elements in
//        the original Matrix.
//
//    ConstSymMatrixView<T> SymMatrixViewOf(const T* m, size_t size, int nlo,
//        uplo, stor)
//    ConstSymMatrixView<T> HermMatrixViewOf(const T* m, size_t size, int nlo,
//        uplo, stor)
//    SymMatrixView<T> SymMatrixViewOf(T* m, size_t size, int nlo, uplo, stor)
//    SymMatrixView<T> HermMatrixViewOf(T* m, size_t size, int nlo, uplo, stor)
//        View the actual memory pointed to by m as a 
//        SymBandMatrix/HermBandMatrix with the given size and storage.
//
// Access Functions
//
//    size_t colsize() const
//    size_t rowsize() const
//    size_t size() const
//    int nlo() const
//    int nhi() const
//        Return the dimensions of the SymMatrix.  
//        (colsize()==rowsize()==size() and nlo()==nhi().)
//
//    T& operator()(int i, int j)
//    T operator()(int i, int j) const
//        Return the (i,j) element of the SymMatrix
//
//    Vector& row(int i, int j1, int j2)
//        Return a portion of the ith row 
//        This range must be a valid range for the requested row
//        and be entirely in the upper or lower band.
//
//    Vector& col(int j, int i1, int i2)
//        Return a portion of the jth column
//        This range must be a valid range for the requested column
//        and be entirely in the upper or lower band.
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
//    void Swap(SymBandMatrix& m1, SymBandMatrix& m2)
//        The SymBandMatrices must be the same size and Hermitianity.
//
// Views of a SymBandMatrix:
//
//    SubMatrix(int i1, int i2, int j1, int j2, int istep=1, int jstep=1)
//    SubBandMatrix(int i1, int i2, int j1, int j2, int nlo, int nhi, 
//        int istep=1, int jstep=1)
//    SubVector(int i, int j, int istep, int jstep, int size)
//        Just like for a regular band matrix except that the 
//        sub-matrix/bandmatrix/vector must be completely contained within 
//        either the upper or lower band.
//
//    SubSymBandMatrix(int i1, int i2, int nlo, int istep)
//    SubSymMatrix(int i1, int i2, int istep)
//        Returns the SymMatrix or SymBandMatrix which runs from i1 to i2 
//        along the diagonal (not including i2) with an optional step, 
//        and includes the off diagonals in the same rows/cols.
//        The first version allows for having fewer off-diagonals than the 
//        original.
//        The parameter nlo may be omitted in which case it is taken to be the 
//        same number of off-diagonals as the original matrix.
//        The second version must have i2-i1 <= nlo+1 (ie. the SymMatrix is 
//        fully contained within the band structure).
//
//    SymDiags(int nlo)
//        Returns a SymBandMatrix which has the same size as the original
//        but fewer off-diagonals.
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
//          sB/hB size nlo
//          ( m(0,0) )
//          ( m(1,0) m(1,1) )
//          ...
//          ( m(nlo,0) ... m(nlo,nlo) )
//          ( m(size-1,size-nlo-1) ... m(size-1,size-1) )
//
//    is >> m
//        Reads m from istream is in the compact format
//        m must already be the correct size for this to work.
//
//    is >> mptr
//        If you do not know the size of the SymBandMatrix to be read, you can
//        use this form where mptr is an auto_ptr to an undefined SymBandMatrix.
//
//
// Division Control Functions:
//
//    m.DivideUsing(dt)
//    where dt is LU, CH, or SV
//     
//    LUD(), CHD(), SVD(), and SymSVD() return the corresponding Divider 
//        classes.
//
//    As for SymMatrix, LU actually does an LDLT decomposition.
//        ie. L(Lower Triangle) * Diagonal * L.Transpose()
//        For non-positive definite matrices, D may have 2x2 blocks along
//        the diagonal, which can then be further decomposed via normal LU
//        decomposition.
//        
//    Remember that Cholskey decomposition is only appropriate if you know 
//        that your HermMatrix is positive definite.  In this case, the 
//        HermmMatrix can be decomposed into L*L.Adjoint().  
//
//    Finally, a difference for the SV Decomposition for Hermitian matrices
//        is that the decomposition can be done with U = Vt.  In this case,
//        the "singular values" are really the eigenvalues of A.
//        The only caveat is that they may be negative, whereas the usual
//        definition of singular values is that they are all positive.
//        Requiring positivity would destroy U = Vt, so it seemed more
//        useful to leave them as the actual eigenvalues.  Just keep that
//        in mind if you use the singular values for anything that expects
//        them to be positive.


#ifndef TMV_SymBandMatrix_H
#define TMV_SymBandMatrix_H

#include "tmv/TMV_BaseSymBandMatrix.h"
#include "tmv/TMV_SymBandMatrixArithFunc.h"
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_DiagMatrix.h"
#include <vector>

namespace tmv {

  template <class T> 
  class GenSymBandMatrix : 
    virtual public AssignableToSymBandMatrix<T>,
    virtual public AssignableToDiagMatrix<T>,
    public BaseMatrix<T>,
    private DivHelper<T>
  {

  public:

    //
    // Constructors
    //

    inline GenSymBandMatrix() {}
    inline GenSymBandMatrix(const GenSymBandMatrix<T>& rhs) {}
    virtual inline ~GenSymBandMatrix() {}

    //
    // Access Functions
    //

    using AssignableToSymMatrix<T>::size;
    inline size_t colsize() const { return size(); }
    inline size_t rowsize() const { return size(); }
    using AssignableToSymMatrix<T>::sym;
    using AssignableToBandMatrix<T>::nlo;
    inline int nhi() const { return nlo(); }

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
      TMVAssert(j1==j2 || okij(i,j1));
      TMVAssert(j1==j2 || okij(i,j2-1));
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
      TMVAssert(i1==i2 || okij(i1,j));
      TMVAssert(i1==i2 || okij(i2-1,j));
      TMVAssert( j-i1<=0 || i2-j<=1 );
      if ((i2-j<=1 && uplo()==Upper) || (j-i1<=0 && uplo()==Lower))
        return ConstVectorView<T>(cptr()+i1*stepi()+j*stepj(),i2-i1,stepi(),
            ct()); 
      else 
        return ConstVectorView<T>(cptr()+i1*stepj()+j*stepi(),i2-i1,stepj(),
            issym()?ct():ConjOf(T,ct()));
    }

    inline ConstVectorView<T> diag() const
    { return ConstVectorView<T>(cptr(),size(),diagstep(),ct()); }

    inline ConstVectorView<T> diag(int i) const
    {
      TMVAssert(i>=-nlo() && i<=nlo());
      TMVAssert(i>=-int(size()) && i<=int(size())); 
      if (i>=0)
        if (uplo()==Upper)
          return ConstVectorView<T>(cptr()+i*stepj(),size()-i,
              diagstep(),ct()); 
        else
          return ConstVectorView<T>(cptr()+i*stepi(),size()-i,
              diagstep(),issym()?ct():ConjOf(T,ct()));
      else
        if (uplo()==Upper)
          return ConstVectorView<T>(cptr()-i*stepj(),size()+i,
              diagstep(),issym()?ct():ConjOf(T,ct()));
        else
          return ConstVectorView<T>(cptr()-i*stepi(),size()+i,
              diagstep(),ct()); 
    }

    inline ConstVectorView<T> diag(int i, int j1, int j2) const
    {
      TMVAssert(i>=-nlo() && i<=nlo());
      TMVAssert(i>=-int(size()) && i<=int(size())); 
      TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size())-std::abs(i));
      if (i>=0)
        if (uplo()==Upper)
          return ConstVectorView<T>(cptr()+i*stepj()+j1*diagstep(),
              j2-j1,diagstep(),ct()); 
        else
          return ConstVectorView<T>(cptr()+i*stepi()+j1*diagstep(),
              j2-j1,diagstep(),issym()?ct():ConjOf(T,ct()));
      else
        if (uplo()==Upper)
          return ConstVectorView<T>(cptr()-i*stepj()+j1*diagstep(),
              j2-j1,diagstep(),issym()?ct():ConjOf(T,ct()));
        else
          return ConstVectorView<T>(cptr()-i*stepi()+j1*diagstep(),
              j2-j1,diagstep(),ct()); 
    }

    template <class T2> 
    inline bool SameAs(const BaseMatrix<T2>& ) const
    { return false; }

    inline bool SameAs(const GenSymBandMatrix<T>& m2) const
    { 
      if (this == &m2) return true;
      else if (cptr()==m2.cptr() && size()==m2.size() && 
          nlo()==m2.nlo() && sym() == m2.sym()) {
        if (uplo() == m2.uplo())
          return (stepi() == m2.stepi() && stepj() == m2.stepj() && 
              ct() == m2.ct());
        else
          return (stepi() == m2.stepj() && stepj() == m2.stepi() && 
              issym() == (ct()==m2.ct()));
      } else return false;
    }

    inline void AssignToM(const MatrixView<RealType(T)>& m2) const
    {
      TMVAssert(IsReal(T()));
      TMVAssert(m2.colsize() == size());
      TMVAssert(m2.rowsize() == size());
      AssignToB(BandMatrixViewOf(m2,nlo(),nlo()));
      if (int(size()) > nlo()+1) {
        m2.UpperTri().OffDiag(nlo()+1).Zero();
        m2.LowerTri().OffDiag(nlo()+1).Zero();
      }
    }

    inline void AssignToM(const MatrixView<ComplexType(T)>& m2) const
    {
      TMVAssert(m2.colsize() == size());
      TMVAssert(m2.rowsize() == size());
      AssignToB(BandMatrixViewOf(m2,nlo(),nlo()));
      if (int(size()) > nlo()+1) {
        m2.UpperTri().OffDiag(nlo()+1).Zero();
        m2.LowerTri().OffDiag(nlo()+1).Zero();
      }
    }

    inline void AssignToD(const DiagMatrixView<RealType(T)>& m2) const
    {
      TMVAssert(IsReal(T()));
      TMVAssert(m2.size() == size());
      TMVAssert(nlo() == 0);
      m2.diag() = diag();
    }

    inline void AssignToD(const DiagMatrixView<ComplexType(T)>& m2) const
    {
      TMVAssert(m2.size() == size());
      TMVAssert(nlo() == 0);
      m2.diag() = diag();
      if (!issym()) m2.diag().Imag().Zero();
    }

    inline void AssignToB(const BandMatrixView<RealType(T)>& m2) const
    {
      TMVAssert(IsReal(T()));
      TMVAssert(m2.colsize() == size());
      TMVAssert(m2.rowsize() == size());
      TMVAssert(m2.nlo() >= nlo());
      TMVAssert(m2.nhi() >= nlo());
      AssignTosB(SymBandMatrixViewOf(m2,Upper,nlo()));
      if (nlo() > 0) m2.Diags(-nlo(),0) = m2.Diags(1,nlo()+1).Transpose();
      if (m2.nlo() > nlo()) m2.Diags(-m2.nlo(),-nlo()).Zero();
      if (m2.nhi() > nlo()) m2.Diags(nlo()+1,m2.nhi()+1).Zero();
    }

    inline void AssignToB(const BandMatrixView<ComplexType(T)>& m2) const
    {
      TMVAssert(m2.colsize() == size());
      TMVAssert(m2.rowsize() == size());
      TMVAssert(m2.nlo() >= nlo());
      TMVAssert(m2.nhi() >= nlo());
      if (issym()) {
        AssignTosB(SymBandMatrixViewOf(m2,Upper,nlo()));
        if (nlo() > 0) m2.Diags(-nlo(),0) = m2.Diags(1,nlo()+1).Transpose();
      } else {
        m2.diag().Imag().Zero();
        AssignTosB(HermBandMatrixViewOf(m2,Upper,nlo()));
        if (nlo() > 0) m2.Diags(-nlo(),0) = m2.Diags(1,nlo()+1).Adjoint();
      }
      if (m2.nlo() > nlo()) m2.Diags(-m2.nlo(),-nlo()).Zero();
      if (m2.nhi() > nlo()) m2.Diags(nlo()+1,m2.nhi()+1).Zero();
    }

    inline void AssignToS(const SymMatrixView<RealType(T)>& m2) const
    {
      TMVAssert(IsReal(T()));
      TMVAssert(m2.size() == size());
      AssignTosB(SymBandMatrixViewOf(m2,nlo()));
      if (int(size()) > nlo()+1) m2.UpperTri().OffDiag(nlo()+1).Zero();
    }

    inline void AssignToS(const SymMatrixView<ComplexType(T)>& m2) const
    {
      TMVAssert(m2.size() == size());
      TMVAssert(m2.sym() == sym());
      if (issym()) AssignTosB(SymBandMatrixViewOf(m2,nlo()));
      else AssignTosB(HermBandMatrixViewOf(m2,nlo()));
      if (int(size()) > nlo()+1) m2.UpperTri().OffDiag(nlo()+1).Zero();
    }

    inline void AssignTosB(const SymBandMatrixView<RealType(T)>& m2) const
    {
      TMVAssert(IsReal(T()));
      TMVAssert(m2.size() == size());
      TMVAssert(m2.nlo() >= nlo());

      if (!SameAs(m2)) m2.UpperBand() = UpperBand();
      if (m2.nlo() > nlo()) m2.Diags(-m2.nlo(),-nlo()).Zero();
    }

    inline void AssignTosB(const SymBandMatrixView<ComplexType(T)>& m2) const
    {
      TMVAssert(m2.size() == size());
      TMVAssert(m2.nlo() >= nlo());
      TMVAssert(IsReal(T()) || m2.sym() == sym());

      if (!SameAs(m2)) m2.UpperBand() = UpperBand();
      if (m2.nlo() > nlo()) m2.Diags(-m2.nlo(),-nlo()).Zero();
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

    bool OKSubBandMatrix(
        int i1, int i2, int j1, int j2, int newnlo, int newnhi,
        int istep, int jstep) const;

    inline ConstBandMatrixView<T> SubBandMatrix(
        int i1, int i2, int j1, int j2, int newnlo, int newnhi) const
    {
      TMVAssert(OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
      TMVAssert(i1+newnlo-j1<=0 || j1+newnhi-i1<=0);
      if ((i1+newnlo-j1<=0 && uplo()==Upper) || 
          (j1+newnhi-i1<=0 && uplo()==Lower))
        return ConstBandMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
            i2-i1,j2-j1,newnlo,newnhi,stepi(),stepj(),diagstep(),
            stor(),ct());
      else
        return ConstBandMatrixView<T>(cptr()+i1*stepj()+j1*stepi(),
            i2-i1,j2-j1,newnlo,newnhi,stepj(),stepi(),diagstep(),
            TransOf(stor()),issym()?ct():ConjOf(T,ct()));
    }

    inline ConstBandMatrixView<T> SubBandMatrix(
        int i1, int i2, int j1, int j2, int newnlo, int newnhi,
        int istep, int jstep) const
    {
      TMVAssert(OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep));
      const StorageType newstor=
      iscm() ? (istep == 1 ? ColMajor : NoMajor) :
      isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
      TMVAssert(i1+newnlo*istep<=j1 || j1+newnhi*jstep<=i1);
      if ((i1+newnlo*istep<=j1 && uplo()==Upper) || 
          (j1+newnhi*jstep<=i1 && uplo()==Lower)) {
        const int newstepi = stepi()*istep;
        const int newstepj = stepj()*jstep;
        return ConstBandMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
            (i2-i1)/istep,(j2-j1)/jstep,newnlo,newnhi,
            newstepi,newstepj,newstepi+newstepj,newstor,ct());
      } else {
        const int newstepi = stepj()*istep;
        const int newstepj = stepi()*jstep;
        return ConstBandMatrixView<T>(cptr()+i1*stepj()+j1*stepi(),
            (i2-i1)/istep,(j2-j1)/jstep,newnlo,newnhi,
            newstepi,newstepj,newstepi+newstepj,
            TransOf(newstor),issym()?ct():ConjOf(T,ct()));
      }
    }

    bool OKSubVector(int i, int j, int istep, int jstep, int n) const;

    inline ConstVectorView<T> SubVector(int i, int j,
        int istep, int jstep, int n) const
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
      return ConstSymMatrixView<T>(cptr()+i1*diagstep(),i2-i1,
          stepi(),stepj(),sym(),uplo(),stor(),ct());
    }

    inline ConstSymMatrixView<T> SubSymMatrix(int i1, int i2,
        int istep) const
    {
      TMVAssert(OKSubSymMatrix(i1,i2,istep));
      return ConstSymMatrixView<T>(cptr()+i1*diagstep(),
          (i2-i1)/istep,istep*stepi(),istep*stepj(),
          sym(),uplo(),istep==1 ? stor() : NoMajor,ct());
    }

    bool OKSubSymBandMatrix(int i1, int i2, int newnlo, int istep) const;

    inline ConstSymBandMatrixView<T> SubSymBandMatrix(
        int i1, int i2, int newnlo) const
    {
      TMVAssert(OKSubSymBandMatrix(i1,i2,newnlo,1));
      return ConstSymBandMatrixView<T>(cptr()+i1*diagstep(),i2-i1,newnlo,
          stepi(),stepj(),diagstep(),sym(),uplo(),stor(),ct());
    }

    inline ConstSymBandMatrixView<T> SubSymBandMatrix(int i1, int i2) const
    { return SubSymBandMatrix(i1,i2,MIN(nlo(),i2-i1-1)); }

    inline ConstSymMatrixView<T> SubSymBandMatrix(
        int i1, int i2, int newnlo, int istep) const
    {
      TMVAssert(OKSubSymBandMatrix(i1,i2,newnlo,istep));
      return ConstSymBandMatrixView<T>(cptr()+i1*(diagstep()),
          (i2-i1)/istep,istep*stepi(),istep*stepj(),istep*diagstep(),
          sym(),uplo(),istep==1 ? stor() : NoMajor,ct());
    }

    inline ConstSymBandMatrixView<T> SymDiags(int newnlo) const
    {
      TMVAssert(newnlo>=0 && newnlo <= nlo());
      return ConstSymBandMatrixView<T>(cptr(),size(),newnlo,
          stepi(),stepj(),diagstep(),sym(),uplo(),stor(),ct());
    }

    inline ConstBandMatrixView<T> Diags(int k1, int k2) const
    {
      TMVAssert(k1>=-nlo() && k1<k2 && k2<=nlo()+1);
      TMVAssert(k2 <= 1 || k1 >= 0);

      if (k1 >= 0) {
        const int newsize = int(size())-k1;
        const int newnhi = k2-k1-1;
        if (uplo() == Upper)
          return ConstBandMatrixView<T>(cptr()+k1*stepj(),
              newsize, newsize, 0, newnhi, stepi(), stepj(), diagstep(),
              stor(), ct());
        else
          return ConstBandMatrixView<T>(cptr()+k1*stepi(),
              newsize, newsize, 0, newnhi, stepj(), stepi(), diagstep(),
              TransOf(stor()), issym()?ct():ConjOf(T,ct()));
      } else {
        const int newsize = int(size())+k2-1;
        const int newnlo = k2-k1-1;
        if (uplo() == Lower)
          return ConstBandMatrixView<T>(cptr()-k2*stepi(),
              newsize, newsize, newnlo, 0, stepi(), stepj(), diagstep(),
              stor(), ct());
        else
          return ConstBandMatrixView<T>(cptr()-k2*stepj(),
              newsize, newsize, newnlo, 0, stepj(), stepi(), diagstep(),
              TransOf(stor()), issym()?ct():ConjOf(T,ct()));
      }
    }

    inline ConstBandMatrixView<T> UpperBand() const
    {
      if (uplo() == Upper)
        return ConstBandMatrixView<T>(cptr(),size(),size(),0,nlo(),
            stepi(),stepj(),diagstep(),stor(),ct());
      else
        return ConstBandMatrixView<T>(cptr(),size(),size(),0,nlo(),
            stepj(),stepi(),diagstep(),TransOf(stor()),
            issym()?ct():ConjOf(T,ct()));
    }

    inline ConstBandMatrixView<T> UpperBandOff() const
    {
      TMVAssert(size()>0);
      TMVAssert(nlo()>0);
      if (uplo() == Upper)
        return ConstBandMatrixView<T>(cptr()+stepj(),size()-1,size()-1,
            0,nlo()-1,stepi(),stepj(),diagstep(),stor(),ct());
      else
        return ConstBandMatrixView<T>(cptr()+stepi(),size()-1,size()-1,
            0,nlo()-1,stepj(),stepi(),diagstep(),TransOf(stor()),
            issym()?ct():ConjOf(T,ct()));
    }

    inline ConstBandMatrixView<T> LowerBand() const
    {
      if (uplo() == Lower)
        return ConstBandMatrixView<T>(cptr(),size(),size(),nlo(),0,
            stepi(),stepj(),diagstep(),stor(),ct());
      else
        return ConstBandMatrixView<T>(cptr(),size(),size(),nlo(),0,
            stepj(),stepi(),diagstep(),TransOf(stor()),
            issym()?ct():ConjOf(T,ct()));
    }

    inline ConstBandMatrixView<T> LowerBandOff() const
    {
      TMVAssert(size()>0);
      TMVAssert(nlo()>0);
      if (uplo() == Lower)
        return ConstBandMatrixView<T>(cptr()+stepi(),size()-1,size()-1,
            nlo()-1,0,stepi(),stepj(),diagstep(),stor(),ct());
      else
        return ConstBandMatrixView<T>(cptr()+stepj(),size()-1,size()-1,
            nlo()-1,0,stepj(),stepi(),diagstep(),TransOf(stor()),
            issym()?ct():ConjOf(T,ct()));
    }

    inline ConstSymBandMatrixView<RealType(T)> Real() const
    {
      return ConstSymBandMatrixView<RealType(T)>(
          reinterpret_cast<const RealType(T)*>(cptr()),size(),nlo(),
          IsReal(T()) ? stepi() : 2*stepi(),
          IsReal(T()) ? stepj() : 2*stepj(),
          IsReal(T()) ? diagstep() : 2*diagstep(),
          Sym, uplo(), IsReal(T()) ? stor() : NoMajor,NonConj);
    }

    inline ConstSymBandMatrixView<RealType(T)> Imag() const
    {
      TMVAssert(IsComplex(T()));
      TMVAssert(issym());
      // The imaginary part of a Hermitian matrix is anti-symmetric
      return ConstSymBandMatrixView<RealType(T)>(
          reinterpret_cast<const RealType(T)*>(cptr())+1,
          size(),nlo(),2*stepi(),2*stepj(),2*diagstep(),
          Sym,uplo(),NoMajor,NonConj);
    }

    // 
    // Views
    //

    inline ConstSymBandMatrixView<T> View() const
    { 
      return ConstSymBandMatrixView<T>(cptr(),size(),nlo(),
          stepi(),stepj(),diagstep(),sym(),uplo(),stor(),ct());
    }

    inline ConstSymBandMatrixView<T> Transpose() const
    { 
      return ConstSymBandMatrixView<T>(cptr(),size(),nlo(),
          stepj(),stepi(),diagstep(),sym(),
          UTransOf(uplo()),TransOf(stor()),ct());
    }

    inline ConstSymBandMatrixView<T> Conjugate() const
    { 
      return ConstSymBandMatrixView<T>(cptr(),size(),nlo(),
          stepi(),stepj(),diagstep(),sym(),uplo(),stor(),ConjOf(T,ct()));
    }

    inline ConstSymBandMatrixView<T> Adjoint() const
    { 
      return ConstSymBandMatrixView<T>(cptr(),size(),nlo(),
          stepj(),stepi(),diagstep(),sym(),
          UTransOf(uplo()),TransOf(stor()),ConjOf(T,ct()));
    }

    inline SymBandMatrixView<T> NonConst() const
    {
      return SymBandMatrixView<T>(const_cast<T*>(cptr()),size(),nlo(),
          stepi(),stepj(),diagstep(),sym(),uplo(),stor(),ct()
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
      TMV_Warning("Calling SymBandMatrix::Norm2 without previously calling DivideUsing(SV)");
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
      TMV_Warning("Calling SymBandMatrix::Condition without previously calling DivideUsing(SV)");
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

    QuotXsB<T,T> QInverse() const;
    inline QuotXsB<T,T> Inverse() const
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

    inline const BandLUDiv<T>& LUD() const
    {
      DivideUsing(LU);
      SetDiv();
      TMVAssert(GetDiv());
      TMVAssert(dynamic_cast<const BandLUDiv<T>*>(GetDiv()));
      return *dynamic_cast<const BandLUDiv<T>*>(GetDiv());
    }

    inline const HermBandCHDiv<T>& CHD() const
    {
      TMVAssert(isherm());
      DivideUsing(CH);
      SetDiv();
      TMVAssert(GetDiv());
      TMVAssert(dynamic_cast<const HermBandCHDiv<T>*>(GetDiv()));
      return *dynamic_cast<const HermBandCHDiv<T>*>(GetDiv());
    }

    inline const HermBandSVDiv<T>& SVD() const
    {
      TMVAssert(isherm());
      DivideUsing(SV);
      SetDiv();
      TMVAssert(GetDiv());
      TMVAssert(dynamic_cast<const HermBandSVDiv<T>*>(GetDiv()));
      return *dynamic_cast<const HermBandSVDiv<T>*>(GetDiv());
    }

    inline const SymBandSVDiv<T>& SymSVD() const
    {
      TMVAssert(IsComplex(T()));
      TMVAssert(issym());
      DivideUsing(SV);
      SetDiv();
      TMVAssert(GetDiv());
      TMVAssert(dynamic_cast<const SymBandSVDiv<T>*>(GetDiv()));
      return *dynamic_cast<const SymBandSVDiv<T>*>(GetDiv());
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
    inline void LDiv(const GenVector<T1>& v1, const VectorView<T0>& v0) const
    { DivHelper<T>::LDiv(v1,v0); }
    template <class T1, class T0> 
    inline void LDiv(const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
    { DivHelper<T>::LDiv(m1,m0); }
    template <class T1, class T0> 
    inline void RDiv(const GenVector<T1>& v1, const VectorView<T0>& v0) const
    { DivHelper<T>::RDiv(v1,v0); }
    template <class T1, class T0> 
    inline void RDiv(const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
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
    virtual int diagstep() const = 0;
    virtual UpLoType uplo() const = 0;
    virtual StorageType stor() const = 0;
    virtual ConjType ct() const = 0;
    virtual inline bool isrm() const { return stor() == RowMajor; }
    virtual inline bool iscm() const { return stor() == ColMajor; }
    virtual inline bool isdm() const { return stor() == DiagMajor; }
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

    inline bool okij(int i, int j) const
    { return (j+nlo() >= i && i+nlo() >= j); }

    void NewDivider() const;
    inline const BaseMatrix<T>& GetMatrix() const { return *this; }

  private :

    inline GenSymBandMatrix<T>& operator=(const GenSymBandMatrix<T>&) 
    { TMVAssert(FALSE); return *this; }

  }; // GenSymBandMatrix

  template <class T, IndexStyle I> 
  class ConstSymBandMatrixView : 
    public GenSymBandMatrix<T>
  {
  public :

    inline ConstSymBandMatrixView(const ConstSymBandMatrixView<T,I>& rhs) :
      itsm(rhs.itsm), itss(rhs.itss), itslo(rhs.itslo),
      itssi(rhs.itssi), itssj(rhs.itssj), itssd(rhs.itssd),
      itssym(rhs.itssym), itsuplo(rhs.itsuplo), itsstor(rhs.itsstor),
      itsct(rhs.itsct) 
    {
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    inline ConstSymBandMatrixView(const GenSymBandMatrix<T>& rhs) :
      itsm(rhs.cptr()), itss(rhs.size()), itslo(rhs.nlo()), 
      itssi(rhs.stepi()), itssj(rhs.stepj()), itssd(rhs.diagstep()),
      itssym(rhs.sym()), itsuplo(rhs.uplo()), itsstor(rhs.stor()),
      itsct(rhs.ct())
    {
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    inline ConstSymBandMatrixView(
        const T* _m, size_t _s, int _lo, int _si, int _sj, int _sd,
        SymType _sym, UpLoType _uplo, StorageType _stor,
        ConjType _ct) : 
      itsm(_m), itss(_s), itslo(_lo), itssi(_si), itssj(_sj), itssd(_sd),
      itssym(_sym), itsuplo(_uplo), itsstor(_stor), itsct(_ct)
    { 
      TMVAssert(_stor==RowMajor ? _sj==1 : _stor==ColMajor ?
          _si==1 : _stor==DiagMajor ? _sd==1 : true); 
      TMVAssert(size()==0 || nlo() < int(size())); 
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    virtual inline ~ConstSymBandMatrixView()
    {
#ifdef TMVDEBUG
      const_cast<const T*&>(itsm) = 0;
#endif
    }

    inline size_t size() const { return itss; }
    inline int nlo() const { return itslo; }
    inline const T* cptr() const { return itsm; }
    inline int stepi() const { return itssi; }
    inline int stepj() const { return itssj; }
    inline int diagstep() const { return itssd; }
    inline SymType sym() const { return itssym; }
    inline UpLoType uplo() const { return itsuplo; }
    inline StorageType stor() const { return itsstor; }
    inline ConjType ct() const { return itsct; }

  protected :

    const T*const itsm;
    const size_t itss;
    const int itslo;
    const int itssi;
    const int itssj;
    const int itssd;

    const SymType itssym;
    const UpLoType itsuplo;
    const StorageType itsstor;
    const ConjType itsct;

  private :

    inline ConstSymBandMatrixView<T,I>& operator=(
        const ConstSymBandMatrixView<T,I>&) 
    { TMVAssert(FALSE); return *this; }

  }; // ConstSymBandMatrixView

  template <class T> 
  class ConstSymBandMatrixView<T,FortranStyle> : 
  public ConstSymBandMatrixView<T,CStyle>
  {
  public :

    inline ConstSymBandMatrixView(
        const ConstSymBandMatrixView<T,FortranStyle>& rhs) :
      ConstSymBandMatrixView<T,CStyle>(rhs) {}

    inline ConstSymBandMatrixView(
        const ConstSymBandMatrixView<T,CStyle>& rhs) :
      ConstSymBandMatrixView<T,CStyle>(rhs) {}

    inline ConstSymBandMatrixView(const GenSymBandMatrix<T>& rhs) :
      ConstSymBandMatrixView<T,CStyle>(rhs) {}

    inline ConstSymBandMatrixView(const T* _m, size_t _s, int _lo,
        int _si, int _sj, int _sd,
        SymType _sym, UpLoType _uplo, StorageType _stor, ConjType _ct) : 
      ConstSymBandMatrixView<T,CStyle>(_m,_s,_lo,_si,_sj,_sd,
          _sym,_uplo,_stor,_ct) {}

    virtual inline ~ConstSymBandMatrixView() {}

    //
    // Access Functions
    //

    inline T operator()(int i, int j) const
    {
      TMVAssert(i>0 && i<=int(this->size()));
      TMVAssert(j>0 && j<=int(this->size()));
      return GenSymBandMatrix<T>::operator()(i-1,j-1);
    }

    inline ConstVectorView<T,FortranStyle> row(int i, int j1, int j2) const 
    { 
      TMVAssert(i>0 && i<=int(this->size()));
      TMVAssert(j1>0 && j1-j2<=0 && j2<=int(this->size()));
      TMVAssert(this->okij(i-1,j1-1));
      TMVAssert(this->okij(i-1,j2-1));
      return GenSymBandMatrix<T>::row(i-1,j1-1,j2);
    }

    inline ConstVectorView<T,FortranStyle> col(int j, int i1, int i2) const
    {
      TMVAssert(j>0 && j<=int(this->size()));
      TMVAssert(i1>0 && i1-i2<=0 && i2<=int(this->size()));
      TMVAssert(this->okij(i1-1,j-1));
      TMVAssert(this->okij(i2-1,j-1));
      return GenSymBandMatrix<T>::col(j-1,i1-1,i2);
    }

    inline ConstVectorView<T,FortranStyle> diag() const
    { return GenSymBandMatrix<T>::diag(); }

    inline ConstVectorView<T,FortranStyle> diag(int i) const
    { return GenSymBandMatrix<T>::diag(i); }

    inline ConstVectorView<T,FortranStyle> diag(int i, int j1, int j2) const
    {
      TMVAssert(i>=-this->nlo() && i<=this->nlo());
      TMVAssert(i>=-int(this->size()) && i<=int(this->size())); 
      TMVAssert(j1>0 && j1-j2<=0 && j2<=int(this->size())-std::abs(i));
      return GenSymBandMatrix<T>::diag(i,j1-1,j2);
    }

    //
    // SubMatrix
    //

    bool OKSubMatrix(int i1, int i2, int j1, int j2, 
        int istep, int jstep) const;

    bool OKSubBandMatrix(int i1, int i2, int j1, int j2, 
        int newnlo, int newnhi, int istep, int jstep) const;

    bool OKSubVector(int i, int j, int istep, int jstep, int n) const;

    bool OKSubSymMatrix(int i1, int i2, int istep) const;

    bool OKSubSymBandMatrix(int i1, int i2, int newnlo, int istep) const;

    inline ConstMatrixView<T,FortranStyle> SubMatrix(
        int i1, int i2, int j1, int j2) const
    {
      TMVAssert(OKSubMatrix(i1,i2,j1,j2,1,1));
      return GenSymBandMatrix<T>::SubMatrix(i1-1,i2,j1-1,j2);
    }

    inline ConstMatrixView<T,FortranStyle> SubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) const
    {
      TMVAssert(OKSubMatrix(i1,i2,j1,j2,istep,jstep));
      return GenSymBandMatrix<T>::SubMatrix(i1-1,i2-1+istep,
          j1-1,j2-1+jstep,istep,jstep);
    }

    inline ConstBandMatrixView<T,FortranStyle> SubBandMatrix(
        int i1, int i2, int j1, int j2, int newnlo, int newnhi) const
    {
      TMVAssert(OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
      return GenSymBandMatrix<T>::SubBandMatrix(i1-1,i2,j1-1,j2,
          newnlo,newnhi);
    }

    inline ConstBandMatrixView<T,FortranStyle> SubBandMatrix(
        int i1, int i2, int j1, int j2, int newnlo, int newnhi,
        int istep, int jstep) const
    {
      TMVAssert(OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep));
      return GenSymBandMatrix<T>::SubBandMatrix(i1-1,i2-1+istep,
          j1-1,j2-1+jstep,newnlo,newnhi,istep,jstep);
    }

    inline ConstVectorView<T,FortranStyle> SubVector(
        int i, int j, int istep, int jstep, int n) const
    {
      TMVAssert(OKSubVector(i,j,istep,jstep,n));
      return GenSymBandMatrix<T>::SubVector(i-1,j-1,istep,jstep,n);
    }

    inline ConstSymMatrixView<T,FortranStyle> SubSymMatrix(
        int i1, int i2) const
    {
      TMVAssert(OKSubSymMatrix(i1,i2,1));
      return GenSymBandMatrix<T>::SubSymMatrix(i1-1,i2);
    }

    inline ConstSymMatrixView<T,FortranStyle> SubSymMatrix(
        int i1, int i2, int istep) const
    {
      TMVAssert(OKSubSymMatrix(i1,i2,istep));
      return GenSymBandMatrix<T>::SubSymMatrix(i1-1,i2-1+istep,istep);
    }

    inline ConstSymBandMatrixView<T,FortranStyle> SubSymBandMatrix(
        int i1, int i2, int newnlo) const
    {
      TMVAssert(OKSubSymBandMatrix(i1,i2,newnlo,1));
      return GenSymBandMatrix<T>::SubSymBandMatrix(i1-1,i2,newnlo);
    }

    inline ConstSymBandMatrixView<T,FortranStyle> SubSymBandMatrix(
        int i1, int i2) const
    { return SubSymBandMatrix(i1,i2,MIN(this->nlo(),i2-i1)); }

    inline ConstSymBandMatrixView<T,FortranStyle> SubSymBandMatrix(
        int i1, int i2, int newnlo, int istep) const
    {
      TMVAssert(OKSubSymBandMatrix(i1,i2,newnlo,istep));
      return GenSymBandMatrix<T>::SubSymBandMatrix(i1-1,i2-1+istep,
          newnlo,istep);
    }

    inline ConstSymBandMatrixView<T,FortranStyle> SymDiags(int newnlo) const
    { return GenSymBandMatrix<T>::SymDiags(newnlo); }

    inline ConstBandMatrixView<T,FortranStyle> Diags(int k1, int k2) const
    { return GenSymBandMatrix<T>::Diags(k1,k2); }

    inline ConstBandMatrixView<T,FortranStyle> UpperBand() const
    { return GenSymBandMatrix<T>::UpperBand(); }

    inline ConstBandMatrixView<T,FortranStyle> UpperBandOff() const
    { return GenSymBandMatrix<T>::UpperBandOff(); }

    inline ConstBandMatrixView<T,FortranStyle> LowerBand() const
    { return GenSymBandMatrix<T>::LowerBand(); }

    inline ConstBandMatrixView<T,FortranStyle> LowerBandOff() const
    { return GenSymBandMatrix<T>::LowerBandOff(); }

    inline ConstSymBandMatrixView<RealType(T),FortranStyle> Real() const
    { return GenSymBandMatrix<T>::Real(); }

    inline ConstSymBandMatrixView<RealType(T),FortranStyle> Imag() const
    { return GenSymBandMatrix<T>::Imag(); }

    // 
    // Views
    //

    inline ConstSymBandMatrixView<T,FortranStyle> View() const
    { return GenSymBandMatrix<T>::View(); }

    inline ConstSymBandMatrixView<T,FortranStyle> Transpose() const
    { return GenSymBandMatrix<T>::Transpose(); }

    inline ConstSymBandMatrixView<T,FortranStyle> Conjugate() const
    { return GenSymBandMatrix<T>::Conjugate(); }

    inline ConstSymBandMatrixView<T,FortranStyle> Adjoint() const
    { return GenSymBandMatrix<T>::Adjoint(); }

    inline SymBandMatrixView<T,FortranStyle> NonConst() const
    { return GenSymBandMatrix<T>::NonConst(); }

  private :

    inline ConstSymBandMatrixView<T,FortranStyle>& operator=(
        const ConstSymBandMatrixView<T,FortranStyle>&) 
    { TMVAssert(FALSE); return *this; }

  }; // FortranStyle ConstSymBandMatrixView

  template <class T, IndexStyle I> 
  class SymBandMatrixView : 
    public GenSymBandMatrix<T>
  {

  public:

    //
    // Constructors
    //

    inline SymBandMatrixView(const SymBandMatrixView<T,I>& rhs) : 
      itsm(rhs.itsm), itss(rhs.itss), itslo(rhs.itslo), 
      itssi(rhs.itssi), itssj(rhs.itssj), itssd(rhs.itssd),
      itssym(rhs.itssym), itsuplo(rhs.itsuplo), itsstor(rhs.itsstor),
      itsct(rhs.itsct) DEFFIRSTLAST(rhs.first,rhs.last) 
    {
#ifdef XTEST
      //TMVAssert(this->HermOK());
#endif
    }

    inline SymBandMatrixView(
        T* _m, size_t _s, int _lo, int _si, int _sj, int _sd,
        SymType _sym, UpLoType _uplo, StorageType _stor, ConjType _ct 
        PARAMFIRSTLAST(T) ) :
      itsm(_m), itss(_s), itslo(_lo), itssi(_si), itssj(_sj), itssd(_sd),
      itssym(_sym), itsuplo(_uplo), itsstor(_stor), itsct(_ct)
                                                              DEFFIRSTLAST(_first,_last)
    {
      TMVAssert(_stor==RowMajor ? _sj==1 : _stor==ColMajor ?
          _si==1 : _stor==DiagMajor ? _sd==1 : true); 
      TMVAssert(size()==0 || nlo() < int(size())); 
#ifdef XTEST
      //TMVAssert(this->HermOK());
#endif
    }

    virtual inline ~SymBandMatrixView()
    {
#ifdef TMVDEBUG
      const_cast<T*&>(itsm) = 0;
#endif
    }

    //
    // Op=
    //

    inline const SymBandMatrixView<T,I>& operator=(
        const SymBandMatrixView<T,I>& m2) const
    { 
      TMVAssert(size() == m2.size());
      TMVAssert(IsReal(T()) || m2.sym() == sym());
      TMVAssert(nlo() >= m2.nlo());
      if (!SameAs(m2)) {
        UpperBand() = m2.UpperBand();
        if (nlo() > m2.nlo()) Diags(-nlo(),-m2.nlo()).Zero();
      }
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    inline const SymBandMatrixView<T,I>& operator=(
        const SymBandMatrixView<T,I>& m2) 
    { 
      TMVAssert(size() == m2.size());
      TMVAssert(IsReal(T()) || m2.sym() == sym());
      TMVAssert(nlo() >= m2.nlo());
      if (!SameAs(m2)) {
        UpperBand() = m2.UpperBand();
        if (nlo() > m2.nlo()) Diags(-nlo(),-m2.nlo()).Zero();
      }
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    inline const SymBandMatrixView<T,I>& operator=(
        const GenSymBandMatrix<RealType(T)>& m2) const
    { 
      TMVAssert(size() == m2.size());
      TMVAssert(nlo() >= m2.nlo());
      if (!SameAs(m2)) m2.AssignTosB(*this); 
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this; 
    }

    inline const SymBandMatrixView<T,I>& operator=(
        const GenSymBandMatrix<ComplexType(T)>& m2) const
    { 
      TMVAssert(IsComplex(T()));
      TMVAssert(size() == m2.size());
      TMVAssert(m2.sym() == sym());
      TMVAssert(nlo() >= m2.nlo());
      if (!SameAs(m2)) m2.AssignTosB(*this); 
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this; 
    }

    template <class T2> 
    inline const SymBandMatrixView<T,I>& operator=(
        const GenSymBandMatrix<T2>& m2) const
    { 
      TMVAssert(size() == m2.size());
      TMVAssert(IsComplex(T()) || IsReal(T2()));
      TMVAssert(IsReal(T2()) || m2.sym() == sym());
      TMVAssert(nlo() >= m2.nlo());
      UpperBand() = m2.UpperBand();
      if (nlo() > m2.nlo()) Diags(-nlo(),-m2.nlo()).Zero();
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this; 
    }

    inline const SymBandMatrixView<T,I>& operator=(T x) const 
    { 
      TMVAssert(issym() || IMAG(x) == RealType(T)(0));
      return SetToIdentity(x); 
    }

    inline const SymBandMatrixView<T,I>& operator=(
        const AssignableToSymBandMatrix<RealType(T)>& m2) const
    { 
      TMVAssert(size() == m2.size());
      m2.AssignTosB(View());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    inline const SymBandMatrixView<T,I>& operator=(
        const AssignableToSymBandMatrix<ComplexType(T)>& m2) const
    { 
      TMVAssert(IsComplex(T()));
      TMVAssert(size() == m2.size());
      TMVAssert(m2.sym() == sym());
      m2.AssignTosB(View());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    inline const SymBandMatrixView<T,I>& operator=(
        const GenDiagMatrix<RealType(T)>& m2) const
    { 
      TMVAssert(size() == m2.size());
      m2.AssignToD(DiagMatrixViewOf(diag()));
      if (nlo() > 0) UpperBandOff().Zero();
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    inline const SymBandMatrixView<T,I>& operator=(
        const GenDiagMatrix<ComplexType(T)>& m2) const
    { 
      TMVAssert(IsComplex(T()));
      TMVAssert(issym());
      TMVAssert(size() == m2.size());
      m2.AssignToD(DiagMatrixViewOf(diag()));
      if (nlo() > 0) UpperBandOff().Zero();
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
      TMVAssert(okij(i,j));
      return ref(i,j); 
    }

    inline VectorView<T> row(int i, int j1, int j2) const 
    { 
      TMVAssert(i>=0 && i<int(size()));
      TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size()));
      TMVAssert(j1==j2 || okij(i,j1));
      TMVAssert(j1==j2 || okij(i,j2-1));
      TMVAssert( i-j1<=0 || j2-i<=1 );
      if ((i-j1<=0 && uplo()==Upper) || (j2-i<=1 && uplo()==Lower))
        return VectorView<T>(ptr()+i*stepi()+j1*stepj(),j2-j1,stepj(),
            ct() FIRSTLAST); 
      else
        return VectorView<T>(ptr()+i*stepj()+j1*stepi(),j2-j1,stepi(),
            issym()?ct():ConjOf(T,ct()) FIRSTLAST);
    }

    inline VectorView<T> col(int j, int i1, int i2) const
    {
      TMVAssert(j>=0 && j<int(size()));
      TMVAssert(i1>=0 && i1-i2<=0 && i2<=int(size()));
      TMVAssert(i1==i2 || okij(i1,j));
      TMVAssert(i1==i2 || okij(i2-1,j));
      TMVAssert( j-i1<=0 || i2-j<=1 );
      if ((i2-j<=1 && uplo()==Upper) || (j-i1<=0 && uplo()==Lower))
        return VectorView<T>(ptr()+i1*stepi()+j*stepj(),i2-i1,stepi(),
            ct() FIRSTLAST); 
      else 
        return VectorView<T>(ptr()+i1*stepj()+j*stepi(),i2-i1,stepj(),
            issym()?ct():ConjOf(T,ct()) FIRSTLAST);
    }

    inline VectorView<T> diag() const
    { return VectorView<T>(ptr(),size(),stepi()+stepj(),ct() FIRSTLAST); }

    inline VectorView<T> diag(int i) const
    {
      TMVAssert(i>=-nlo() && i<=nlo());
      TMVAssert(i>=-int(size()) && i<=int(size())); 
      if (i>=0)
        if (uplo()==Upper)
          return VectorView<T>(ptr()+i*stepj(),size()-i,
              diagstep(),ct() FIRSTLAST); 
        else
          return VectorView<T>(ptr()+i*stepi(),size()-i,
              diagstep(),issym()?ct():ConjOf(T,ct()) FIRSTLAST);
      else
        if (uplo()==Upper)
          return VectorView<T>(ptr()-i*stepj(),size()+i,
              diagstep(),issym()?ct():ConjOf(T,ct()) FIRSTLAST);
        else
          return VectorView<T>(ptr()-i*stepi(),size()+i,
              diagstep(),ct() FIRSTLAST); 
    }

    inline VectorView<T> diag(int i, int j1, int j2) const
    {
      TMVAssert(i>=-nlo() && i<=nlo());
      TMVAssert(i>=-int(size()) && i<=int(size())); 
      TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size())-std::abs(i));
      if (i>=0)
        if (uplo()==Upper)
          return VectorView<T>(ptr()+i*stepj()+j1*diagstep(),j2-j1,
              diagstep(),ct() FIRSTLAST); 
        else
          return VectorView<T>(ptr()+i*stepi()+j1*diagstep(),j2-j1,
              diagstep(),issym()?ct():ConjOf(T,ct()) FIRSTLAST);
      else
        if (uplo()==Upper)
          return VectorView<T>(ptr()-i*stepj()+j1*diagstep(),j2-j1,
              diagstep(),issym()?ct():ConjOf(T,ct()) FIRSTLAST);
        else
          return VectorView<T>(ptr()-i*stepi()+j1*diagstep(),j2-j1,
              diagstep(),ct() FIRSTLAST); 
    }

    //
    // Modifying Functions
    //

    inline const SymBandMatrixView<T,I>& Zero() const 
    { UpperBand().Zero(); return *this; }

    inline const SymBandMatrixView<T,I>& SetAllTo(T x) const
    { 
      TMVAssert(IMAG(x)==RealType(T)(0) || this->issym());
      UpperBand().SetAllTo(x); return *this; 
    }

    inline const SymBandMatrixView<T,I>& Clip(RealType(T) thresh) const
    { UpperBand().Clip(thresh); return *this; }

    inline const SymBandMatrixView<T,I>& ConjugateSelf() const
    { if (IsComplex(T())) UpperBand().ConjugateSelf(); return *this; }

    inline const SymBandMatrixView<T,I>& TransposeSelf() const
    { if (!this->issym()) UpperBand().ConjugateSelf(); return *this; }

    inline const SymBandMatrixView<T,I>& SetToIdentity(T x=T(1)) const
    { 
      TMVAssert(IMAG(x)==RealType(T)(0) || this->issym());
      Zero(); diag().SetAllTo(x); return *this; 
    }

    //
    // SubBandMatrix
    //

    inline MatrixView<T> SubMatrix(int i1, int i2, int j1, int j2) const
    {
      TMVAssert(GenSymBandMatrix<T>::OKSubMatrix(i1,i2,j1,j2,1,1));
      TMVAssert(i2-1<=j1 || j2-1<=i1);
      if ((i2-1<=j1 && uplo()==Upper) || (j2-1<=i1 && uplo()==Lower))
        return MatrixView<T>(ptr()+i1*stepi()+j1*stepj(),
            i2-i1,j2-j1,stepi(),stepj(),stor(),ct() FIRSTLAST);
      else
        return MatrixView<T>(ptr()+i1*stepj()+j1*stepi(),
            i2-i1,j2-j1,stepj(),stepi(),TransOf(stor()),
            this->issym()?ct():ConjOf(T,ct()) FIRSTLAST);
    }

    inline MatrixView<T> SubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) const
    {
      const StorageType newstor =
      iscm() ? (istep == 1 ? ColMajor : NoMajor) :
      isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
      TMVAssert(GenSymBandMatrix<T>::OKSubMatrix(i1,i2,j1,j2,istep,jstep));
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

    inline BandMatrixView<T> SubBandMatrix(
        int i1, int i2, int j1, int j2, int newnlo, int newnhi) const
    {
      TMVAssert(GenSymBandMatrix<T>::OKSubBandMatrix(
            i1,i2,j1,j2,newnlo,newnhi,1,1));
      TMVAssert(i1+newnlo-j1<=0 || j1+newnhi-i1<=0);
      if ((i1+newnlo-j1<=0 && uplo()==Upper) || 
          (j1+newnhi-i1<=0 && uplo()==Lower))
        return BandMatrixView<T>(ptr()+i1*stepi()+j1*stepj(),
            i2-i1,j2-j1,newnlo,newnhi,stepi(),stepj(),diagstep(),
            stor(),ct() FIRSTLAST);
      else
        return BandMatrixView<T>(ptr()+i1*stepj()+j1*stepi(),
            i2-i1,j2-j1,newnlo,newnhi,stepj(),stepi(),diagstep(),
            TransOf(stor()), this->issym()?ct():ConjOf(T,ct()) FIRSTLAST);
    }

    inline MatrixView<T> SubBandMatrix(
        int i1, int i2, int j1, int j2, int newnlo, int newnhi,
        int istep, int jstep) const
    {
      TMVAssert(GenSymBandMatrix<T>::OKSubBandMatrix(
            i1,i2,j1,j2,newnlo,newnhi,istep,jstep));
      const StorageType newstor =
      iscm() ? (istep == 1 ? ColMajor : NoMajor) :
      isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
      TMVAssert(i1+newnlo*istep<=j1 || j1+newnhi*jstep<=i1);
      if ((i1+newnlo*istep<=j1 && uplo()==Upper) || 
          (j1+newnhi*jstep<=i1 && uplo()==Lower)) {
        const int newstepi = stepi()*istep;
        const int newstepj = stepj()*jstep;
        return BandMatrixView<T>(ptr()+i1*stepi()+j1*stepj(),
            (i2-i1)/istep,(j2-j1)/jstep,newnlo,newnhi,
            newstepi,newstepj,newstepi+newstepj,newstor,ct() FIRSTLAST);
      } else {
        const int newstepi = stepj()*istep;
        const int newstepj = stepi()*jstep;
        return BandMatrixView<T>(ptr()+i1*stepj()+j1*stepi(),
            (i2-i1)/istep,(j2-j1)/jstep,newnlo,newnhi,
            newstepi,newstepj,newstepi+newstepj,
            TransOf(newstor),this->issym()?ct():ConjOf(T,ct()) FIRSTLAST);
      }
    }

    inline VectorView<T> SubVector(
        int i, int j, int istep, int jstep, int n) const
    {
      TMVAssert(GenSymBandMatrix<T>::OKSubVector(i,j,istep,jstep,n));
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
      TMVAssert(GenSymBandMatrix<T>::OKSubSymMatrix(i1,i2,1));
      return SymMatrixView<T,I>(ptr()+i1*(stepi()+stepj()),i2-i1,
          stepi(),stepj(),sym(),uplo(),stor(),ct() FIRSTLAST);
    }

    inline SymMatrixView<T,I> SubSymMatrix(int i1, int i2, int istep) const
    {
      TMVAssert(GenSymBandMatrix<T>::OKSubSymMatrix(i1,i2,istep));
      return SymMatrixView<T,I>(ptr()+i1*(stepi()+stepj()),(i2-i1)/istep,
          istep*stepi(),istep*stepj(),sym(),uplo(),
          istep==1 ? stor() : NoMajor,ct() FIRSTLAST);
    }

    inline SymBandMatrixView<T,I> SubSymBandMatrix(
        int i1, int i2, int newnlo) const
    {
      TMVAssert(GenSymBandMatrix<T>::OKSubSymBandMatrix(i1,i2,newnlo,1));
      return SymBandMatrixView<T,I>(ptr()+i1*diagstep(),i2-i1,newnlo,
          stepi(),stepj(),diagstep(),sym(),uplo(),stor(),ct() FIRSTLAST);
    }

    inline SymBandMatrixView<T,I> SubSymBandMatrix(int i1, int i2) const
    { return SubSymBandMatrix(i1,i2,MIN(this->nlo(),i2-i1-1)); }

    inline SymBandMatrixView<T,I> SubSymBandMatrix(
        int i1, int i2, int newnlo, int istep) const
    {
      TMVAssert(GenSymBandMatrix<T>::OKSubSymBandMatrix(i1,i2,newnlo,istep));
      return SymBandMatrixView<T,I>(ptr()+i1*diagstep(),(i2-i1)/istep,newnlo,
          istep*stepi(),istep*stepj(),istep*diagstep(),
          sym(),uplo(),istep==1 ? stor() : NoMajor,ct() FIRSTLAST);
    }

    inline SymBandMatrixView<T,I> SymDiags(int newnlo) const
    {
      TMVAssert(newnlo>=0 && newnlo <= nlo());
      return SymBandMatrixView<T,I>(ptr(),size(),newnlo,
          stepi(),stepj(),diagstep(),sym(),uplo(),stor(),ct() FIRSTLAST);
    }

    inline BandMatrixView<T,I> Diags(int k1, int k2) const
    {
      TMVAssert(k1>=-nlo() && k1<k2 && k2<=nlo()+1);
      TMVAssert(k2 <= 1 || k1 >= 0);

      if (k1 >= 0) {
        const int newsize = int(size())-k1;
        const int newnhi = k2-k1-1;
        if (uplo() == Upper)
          return BandMatrixView<T>(ptr()+k1*stepj(),
              newsize, newsize, 0, newnhi, stepi(), stepj(), diagstep(),
              stor(), ct() FIRSTLAST );
        else
          return BandMatrixView<T>(ptr()+k1*stepi(),
              newsize, newsize, 0, newnhi, stepj(), stepi(), diagstep(),
              TransOf(stor()), issym()?ct():ConjOf(T,ct()) FIRSTLAST );
      } else {
        const int newsize = int(size())+k2-1;
        const int newnlo = k2-k1-1;
        if (uplo() == Lower)
          return BandMatrixView<T>(ptr()-k2*stepi(),
              newsize, newsize, newnlo, 0, stepi(), stepj(), diagstep(),
              stor(), ct() FIRSTLAST );
        else
          return BandMatrixView<T>(ptr()-k2*stepj(),
              newsize, newsize, newnlo, 0, stepj(), stepi(), diagstep(),
              TransOf(stor()), issym()?ct():ConjOf(T,ct()) FIRSTLAST );
      }
    }

    inline BandMatrixView<T,I> UpperBand() const
    {
      if (uplo() == Upper)
        return BandMatrixView<T>(ptr(),size(),size(),0,nlo(),
            stepi(),stepj(),diagstep(),stor(),ct() FIRSTLAST);
      else
        return BandMatrixView<T>(ptr(),size(),size(),0,nlo(),
            stepj(),stepi(),diagstep(),TransOf(stor()),
            this->issym()?ct():ConjOf(T,ct()) FIRSTLAST);
    }

    inline BandMatrixView<T,I> UpperBandOff() const
    {
      TMVAssert(size()>0);
      TMVAssert(nlo()>0);
      if (uplo() == Upper)
        return BandMatrixView<T>(ptr()+stepj(),size()-1,size()-1,
            0,nlo()-1,stepi(),stepj(),diagstep(),stor(),ct() FIRSTLAST);
      else
        return BandMatrixView<T>(ptr()+stepi(),size()-1,size()-1,
            0,nlo()-1,stepj(),stepi(),diagstep(),TransOf(stor()),
            this->issym()?ct():ConjOf(T,ct()) FIRSTLAST);
    }

    inline BandMatrixView<T,I> LowerBand() const
    {
      if (uplo() == Lower)
        return BandMatrixView<T>(ptr(),size(),size(),nlo(),0,
            stepi(),stepj(),diagstep(),stor(),ct() FIRSTLAST);
      else
        return BandMatrixView<T>(ptr(),size(),size(),nlo(),0,
            stepj(),stepi(),diagstep(),TransOf(stor()),
            this->issym()?ct():ConjOf(T,ct()) FIRSTLAST);
    }

    inline BandMatrixView<T,I> LowerBandOff() const
    {
      TMVAssert(size()>0);
      TMVAssert(nlo()>0);
      if (uplo() == Lower)
        return BandMatrixView<T>(ptr()+stepi(),size()-1,size()-1,
            nlo()-1,0,stepi(),stepj(),diagstep(),stor(),ct() FIRSTLAST);
      else
        return BandMatrixView<T>(ptr()+stepj(),size()-1,size()-1,
            nlo()-1,0,stepj(),stepi(),diagstep(),TransOf(stor()),
            this->issym()?ct():ConjOf(T,ct()) FIRSTLAST);
    }

    inline SymBandMatrixView<T,I> View() const
    { return *this; }

    inline SymBandMatrixView<T,I> Transpose() const
    {
      return SymBandMatrixView<T,I>(ptr(),size(),nlo(),
          stepj(),stepi(),diagstep(),
          sym(),UTransOf(uplo()),TransOf(stor()),ct() FIRSTLAST);
    }

    inline SymBandMatrixView<T,I> Conjugate() const
    {
      return SymBandMatrixView<T,I>(ptr(),size(),nlo(),
          stepi(),stepj(),diagstep(),
          sym(),uplo(),stor(),ConjOf(T,ct()) FIRSTLAST);
    }

    inline SymBandMatrixView<T,I> Adjoint() const
    {
      return SymBandMatrixView<T,I>(ptr(),size(),nlo(),
          stepj(),stepi(),diagstep(),
          sym(),UTransOf(uplo()),TransOf(stor()),ConjOf(T,ct()) 
          FIRSTLAST);
    }

    inline SymBandMatrixView<RealType(T),I> Real() const
    {
      return SymBandMatrixView<RealType(T)>(
          reinterpret_cast<RealType(T)*>(ptr()),size(),nlo(),
          IsReal(T()) ? stepi() : 2*stepi(),
          IsReal(T()) ? stepj() : 2*stepj(),
          IsReal(T()) ? diagstep() : 2*diagstep(),
          sym(),uplo(), IsReal(T()) ? stor() : NoMajor, NonConj
#ifdef TMVFLDEBUG
          ,reinterpret_cast<const RealType(T)*>(first)
          ,reinterpret_cast<const RealType(T)*>(last)
#endif
          );
    }

    inline SymBandMatrixView<RealType(T),I> Imag() const
    {
      TMVAssert(IsComplex(T()));
      TMVAssert(this->issym());
      return SymBandMatrixView<RealType(T)>(
          reinterpret_cast<RealType(T)*>(ptr())+1,
          size(),nlo(),2*stepi(),2*stepj(),2*diagstep(),
          sym(),uplo(),NoMajor, NonConj
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
    inline int nlo() const { return itslo; }
    inline const T* cptr() const { return itsm; }
    inline T* ptr() const { return itsm; }
    inline int stepi() const { return itssi; }
    inline int stepj() const { return itssj; }
    inline int diagstep() const { return itssd; }
    inline SymType sym() const { return itssym; }
    inline UpLoType uplo() const { return itsuplo; }
    inline StorageType stor() const { return itsstor; }
    inline ConjType ct() const { return itsct; }
    using GenSymBandMatrix<T>::issym;
    using GenSymBandMatrix<T>::iscm;
    using GenSymBandMatrix<T>::isrm;

    reference ref(int i, int j) const;

  protected :

    T*const itsm;
    const size_t itss;
    const int itslo;
    const int itssi;
    const int itssj;
    const int itssd;

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

    using GenSymBandMatrix<T>::okij;

  }; // SymBandMatrixView

  template <class T> 
  class SymBandMatrixView<T,FortranStyle> : 
  public SymBandMatrixView<T,CStyle>
  {

  public:

    //
    // Constructors
    //

    inline SymBandMatrixView(const SymBandMatrixView<T,FortranStyle>& rhs) : 
      SymBandMatrixView<T,CStyle>(rhs) {}

    inline SymBandMatrixView(const SymBandMatrixView<T,CStyle>& rhs) : 
      SymBandMatrixView<T,CStyle>(rhs) {}

    inline SymBandMatrixView(
        T* _m, size_t _s, int _lo, int _si, int _sj, int _sd,
        SymType _sym, UpLoType _uplo, StorageType _stor, ConjType _ct 
        PARAMFIRSTLAST(T) ) :
      SymBandMatrixView<T,CStyle>(_m,_s,_lo,_si,_sj,_sd,
          _sym,_uplo,_stor,_ct FIRSTLAST1(_first,_last) ) {}

    virtual inline ~SymBandMatrixView() {} 

    //
    // Op=
    //

    inline const SymBandMatrixView<T,FortranStyle>& operator=(
        const SymBandMatrixView<T,FortranStyle>& m2) const
    { SymBandMatrixView<T,CStyle>::operator=(m2); return *this; }

    inline const SymBandMatrixView<T,FortranStyle>& operator=(
        const SymBandMatrixView<T,FortranStyle>& m2) 
    { SymBandMatrixView<T,CStyle>::operator=(m2); return *this; }

    inline const SymBandMatrixView<T,FortranStyle>& operator=(
        const SymBandMatrixView<T,CStyle>& m2) const
    { SymBandMatrixView<T,CStyle>::operator=(m2); return *this; }

    inline const SymBandMatrixView<T,FortranStyle>& operator=(
        const SymBandMatrixView<T,CStyle>& m2)
    { SymBandMatrixView<T,CStyle>::operator=(m2); return *this; }

    inline const SymBandMatrixView<T,FortranStyle>& operator=(
        const GenSymBandMatrix<RealType(T)>& m2) const
    { SymBandMatrixView<T,CStyle>::operator=(m2); return *this; }

    inline const SymBandMatrixView<T,FortranStyle>& operator=(
        const GenSymBandMatrix<ComplexType(T)>& m2) const
    { SymBandMatrixView<T,CStyle>::operator=(m2); return *this; }

    template <class T2> 
    inline const SymBandMatrixView<T,FortranStyle>& operator=(
        const GenSymBandMatrix<T2>& m2) const
    { SymBandMatrixView<T,CStyle>::operator=(m2); return *this; }

    inline const SymBandMatrixView<T,FortranStyle>& operator=(T x) const 
    { SymBandMatrixView<T,CStyle>::operator=(x); return *this; }

    inline const SymBandMatrixView<T,FortranStyle>& operator=(
        const AssignableToSymBandMatrix<RealType(T)>& m2) const
    { SymBandMatrixView<T,CStyle>::operator=(m2); return *this; }

    inline const SymBandMatrixView<T,FortranStyle>& operator=(
        const AssignableToSymBandMatrix<ComplexType(T)>& m2) const
    { SymBandMatrixView<T,CStyle>::operator=(m2); return *this; }

    inline const SymBandMatrixView<T,FortranStyle>& operator=(
        const GenDiagMatrix<RealType(T)>& m2) const
    { SymBandMatrixView<T,CStyle>::operator=(m2); return *this; }

    inline const SymBandMatrixView<T,FortranStyle>& operator=(
        const GenDiagMatrix<ComplexType(T)>& m2) const
    { SymBandMatrixView<T,CStyle>::operator=(m2); return *this; }

    inline const SymBandMatrixView<T,FortranStyle>& operator=(
        const GenBandMatrix<RealType(T)>& m2) const
    { SymBandMatrixView<T,CStyle>::operator=(m2); return *this; }

    inline const SymBandMatrixView<T,FortranStyle>& operator=(
        const GenBandMatrix<ComplexType(T)>& m2) const
    { SymBandMatrixView<T,CStyle>::operator=(m2); return *this; }

    //
    // Access
    //

    inline RefType(T) operator()(int i,int j) const 
    {
      TMVAssert(i>0 && i<=int(this->size()));
      TMVAssert(j>0 && j<=int(this->size()));
      TMVAssert(this->okij(i-1,j-1));
      return SymBandMatrixView<T,CStyle>::ref(i-1,j-1); 
    }

    inline VectorView<T,FortranStyle> row(int i, int j1, int j2) const 
    { 
      TMVAssert(i>0 && i<=int(this->size()));
      TMVAssert(j1>0 && j1-j2<=0 && j2<=int(this->size()));
      TMVAssert(this->okij(i-1,j1-1));
      TMVAssert(this->okij(i-1,j2-1));
      return SymBandMatrixView<T,CStyle>::row(i-1,j1-1,j2);
    }

    inline VectorView<T,FortranStyle> col(int j, int i1, int i2) const
    {
      TMVAssert(j>0 && j<=int(this->size()));
      TMVAssert(i1>0 && i1-i2<=0 && i2<=int(this->size()));
      TMVAssert(this->okij(i1-1,j-1));
      TMVAssert(this->okij(i2-1,j-1));
      return SymBandMatrixView<T,CStyle>::col(j-1,i1-1,i2);
    }

    inline VectorView<T,FortranStyle> diag() const
    { return SymBandMatrixView<T,CStyle>::diag(); }

    inline VectorView<T,FortranStyle> diag(int i) const
    { return SymBandMatrixView<T,CStyle>::diag(i); }

    inline VectorView<T,FortranStyle> diag(int i, int j1, int j2) const
    {
      TMVAssert(i>=-this->nlo() && i<=this->nlo());
      TMVAssert(i>=-int(this->size()) && i<=int(this->size())); 
      TMVAssert(j1>0 && j1-j2<=0 && j2<=int(this->size())-std::abs(i));
      return SymBandMatrixView<T,CStyle>::diag(i,j1-1,j2); 
    }

    //
    // Modifying Functions
    //

    inline const SymBandMatrixView<T,FortranStyle>& Zero() const 
    { SymBandMatrixView<T,CStyle>::Zero(); return *this; }

    inline const SymBandMatrixView<T,FortranStyle>& SetAllTo(T x) const
    { SymBandMatrixView<T,CStyle>::SetAllTo(x); return *this; }

    inline const SymBandMatrixView<T,FortranStyle>& Clip(
        RealType(T) thresh) const
    { SymBandMatrixView<T,CStyle>::Clip(thresh); return *this; }

    inline const SymBandMatrixView<T,FortranStyle>& ConjugateSelf() const
    { SymBandMatrixView<T,CStyle>::ConjugateSelf(); return *this; }

    inline const SymBandMatrixView<T,FortranStyle>& TransposeSelf() const
    { SymBandMatrixView<T,CStyle>::TransposeSelf(); return *this; }

    inline const SymBandMatrixView<T,FortranStyle>& SetToIdentity(
        T x=T(1)) const
    { SymBandMatrixView<T,CStyle>::SetToIdentity(x); return *this; }


    //
    // SubMatrix
    //

    inline bool OKSubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) const
    {
      return ConstSymBandMatrixView<T,FortranStyle>(*this).OKSubMatrix(
          i1,i2,j1,j2,istep,jstep); 
    }

    inline bool OKSubBandMatrix(
        int i1, int i2, int j1, int j2, int lo, int hi, 
        int istep, int jstep) const
    {
      return ConstSymBandMatrixView<T,FortranStyle>(*this).OKSubBandMatrix(
          i1,i2,j1,j2,lo,hi,istep,jstep); 
    }

    inline bool OKSubVector(int i, int j, int istep, int jstep,
        int n) const
    {
      return ConstSymBandMatrixView<T,FortranStyle>(*this).OKSubVector(
          i,j,istep,jstep,n); 
    }

    inline bool OKSubSymMatrix(int i1, int i2, int istep) const
    {
      return ConstSymBandMatrixView<T,FortranStyle>(*this).OKSubSymMatrix(
          i1,i2,istep); 
    }

    inline bool OKSubSymBandMatrix(int i1, int i2, int lo, int istep) const
    {
      return ConstSymBandMatrixView<T,FortranStyle>(*this).OKSubSymMatrix(
          i1,i2,lo,istep); 
    }

    inline MatrixView<T,FortranStyle> SubMatrix(
        int i1, int i2, int j1, int j2) const
    {
      TMVAssert(OKSubMatrix(i1,i2,j1,j2,1,1));
      return SymBandMatrixView<T,CStyle>::SubMatrix(i1-1,i2,j1-1,j2);
    }

    inline MatrixView<T,FortranStyle> SubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) const
    {
      TMVAssert(OKSubMatrix(i1,i2,j1,j2,istep,jstep));
      return SymBandMatrixView<T,CStyle>::SubMatrix(i1-1,i2-1+istep,
          j1-1,j2-1+jstep,istep,jstep);
    }

    inline BandMatrixView<T,FortranStyle> SubBandMatrix(
        int i1, int i2, int j1, int j2, int lo, int hi) const
    {
      TMVAssert(OKSubMatrix(i1,i2,j1,j2,lo,hi,1,1));
      return SymBandMatrixView<T,CStyle>::SubBandMatrix(
          i1-1,i2,j1-1,j2,lo,hi);
    }

    inline BandMatrixView<T,FortranStyle> SubBandMatrix(
        int i1, int i2, int j1, int j2, int lo, int hi, 
        int istep, int jstep) const
    {
      TMVAssert(OKSubMatrix(i1,i2,j1,j2,lo,hi,istep,jstep));
      return SymBandMatrixView<T,CStyle>::SubBandMatrix(i1-1,i2-1+istep,
          j1-1,j2-1+jstep,lo,hi,istep,jstep);
    }

    inline VectorView<T,FortranStyle> SubVector(
        int i, int j, int istep, int jstep, int n) const
    {
      TMVAssert(OKSubVector(i,j,istep,jstep,n));
      return SymBandMatrixView<T,CStyle>::SubVector(i-1,j-1,istep,jstep,n);
    }

    inline SymMatrixView<T,FortranStyle> SubSymMatrix(int i1, int i2) const
    {
      TMVAssert(OKSubSymMatrix(i1,i2,1));
      return SymBandMatrixView<T,CStyle>::SubSymMatrix(i1-1,i2);
    }

    inline SymMatrixView<T,FortranStyle> SubSymMatrix(
        int i1, int i2, int istep) const
    {
      TMVAssert(OKSubSymMatrix(i1,i2,istep));
      return SymBandMatrixView<T,CStyle>::SubSymMatrix(
          i1-1,i2-1+istep,istep);
    }

    inline SymBandMatrixView<T,FortranStyle> SubSymBandMatrix(
        int i1, int i2, int lo) const
    {
      TMVAssert(OKSubSymBandMatrix(i1,i2,lo,1));
      return SymBandMatrixView<T,CStyle>::SubSymBandMatrix(i1-1,i2,lo);
    }

    inline SymBandMatrixView<T,FortranStyle> SubSymBandMatrix(
        int i1, int i2) const
    { return SubSymBandMatrix(i1,i2,MIN(this->nlo(),i2-i1)); }

    inline SymMatrixView<T,FortranStyle> SubSymMatrix(
        int i1, int i2, int lo, int istep) const
    {
      TMVAssert(OKSubSymMatrix(i1,i2,lo,istep));
      return SymBandMatrixView<T,CStyle>::SubSymMatrix(
          i1-1,i2-1+istep,lo,istep);
    }

    inline SymBandMatrixView<T,FortranStyle> SymDiags(int newnlo) const
    { return SymBandMatrixView<T,CStyle>::SymDiags(newnlo); }

    inline BandMatrixView<T,FortranStyle> Diags(int k1, int k2) const
    { return SymBandMatrixView<T,CStyle>::Diags(k1,k2); }

    inline BandMatrixView<T,FortranStyle> UpperBand() const
    { return SymBandMatrixView<T,CStyle>::UpperBand(); }

    inline BandMatrixView<T,FortranStyle> UpperBandOff() const
    { return SymBandMatrixView<T,CStyle>::UpperBandOff(); }

    inline BandMatrixView<T,FortranStyle> LowerBand() const
    { return SymBandMatrixView<T,CStyle>::LowerBand(); }

    inline BandMatrixView<T,FortranStyle> LowerBandOff() const
    { return SymBandMatrixView<T,CStyle>::LowerBandOff(); }

    inline SymBandMatrixView<T,FortranStyle> View() const
    { return SymBandMatrixView<T,CStyle>::View(); }

    inline SymBandMatrixView<T,FortranStyle> Transpose() const
    { return SymBandMatrixView<T,CStyle>::Transpose(); }

    inline SymBandMatrixView<T,FortranStyle> Conjugate() const
    { return SymBandMatrixView<T,CStyle>::Conjugate(); }

    inline SymBandMatrixView<T,FortranStyle> Adjoint() const
    { return SymBandMatrixView<T,CStyle>::Adjoint(); }

    inline SymBandMatrixView<RealType(T),FortranStyle> Real() const
    { return SymBandMatrixView<T,CStyle>::Real(); }

    inline SymBandMatrixView<RealType(T),FortranStyle> Imag() const
    { return SymBandMatrixView<T,CStyle>::Imag(); }

  }; // FortranStyle SymBandMatrixView

#ifdef XTEST
#ifdef TMVDEBUG
#define XTEST_DEBUG
#endif
#endif

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  class SymBandMatrix : 
    public GenSymBandMatrix<T>
  {

  public:

    //
    // Constructors
    //

#define NEW_SIZE(s,lo) \
    linsize(BandStorageLength(S,s,s,lo,0)), \
    itsm1(new T[linsize]), itss(s), itslo(lo), \
    itssi(S==DiagMajor ? -int(s)+1 : S==RowMajor ? lo : 1), \
    itssj(S==DiagMajor ? int(s) : S==RowMajor ? 1 : lo), \
    itssd(S==DiagMajor ? 1 : lo+1), \
    itsm((S==DiagMajor && U==Lower) ? itsm1.get()-lo*itssi : itsm1.get()) \
    DEFFIRSTLAST(itsm1.get(),itsm1.get()+linsize) 

    inline SymBandMatrix(size_t s, int lo) : NEW_SIZE(s,lo) 
    { 
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
#ifdef TMVDEBUG
      SetAllTo(T(888));
#endif
    }

    inline SymBandMatrix(size_t s, int lo, const T x) : NEW_SIZE(s,lo) 
    { 
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      SetAllTo(x);
    }

    inline SymBandMatrix(size_t s, int lo, const T* vv) : NEW_SIZE(s,lo)
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      std::copy(vv,vv+linsize,itsm1.get());
    }

    inline SymBandMatrix(size_t s, int lo, const std::vector<T>& vv) : 
      NEW_SIZE(s,lo)
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      TMVAssert(vv.size() == linsize);
      std::copy(vv.begin(),vv.end(),itsm1.get());
    }

    inline SymBandMatrix(const SymBandMatrix<T,U,S,I>& rhs) : 
      NEW_SIZE(rhs.size(),rhs.nlo())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      std::copy(rhs.start_mem(),rhs.start_mem()+linsize,itsm1.get());
    }

    template <UpLoType U2, StorageType S2, IndexStyle I2> 
    inline SymBandMatrix(const SymBandMatrix<T,U2,S2,I2>& rhs) : 
      NEW_SIZE(rhs.size(),rhs.nlo())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      if ( (U==U2 && S==S2) ||
          (S!=DiagMajor && (S2==RowMajor || S2==ColMajor) && 
           S!=S2 && U!=U2) )
        std::copy(rhs.start_mem(),rhs.start_mem()+linsize,itsm1.get());
      else if (U==Upper) UpperBand() = rhs.UpperBand();
      else LowerBand() = rhs.LowerBand();
    }

    inline SymBandMatrix(const GenSymBandMatrix<RealType(T)>& rhs) :
      NEW_SIZE(rhs.size(),rhs.nlo())
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      if (rhs.issym()) {
        rhs.AssignTosB(View());
      } else {
        if (U==Upper) 
          UpperBand() = rhs.UpperBand();
        else 
          LowerBand() = rhs.LowerBand();
      }
    }

    inline SymBandMatrix(const GenSymBandMatrix<ComplexType(T)>& rhs) :
      NEW_SIZE(rhs.size(),rhs.nlo())
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(IsComplex(T()));
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      if (rhs.issym()) {
        rhs.AssignTosB(View());
      } else {
        if (U==Upper) 
          UpperBand() = rhs.UpperBand();
        else 
          LowerBand() = rhs.LowerBand();
      }
    }

    template <class T2> 
    inline SymBandMatrix(const GenSymBandMatrix<T2>& rhs) :
      NEW_SIZE(rhs.size(),rhs.nlo())
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      if (U==Upper) UpperBand() = rhs.UpperBand();
      else LowerBand() = rhs.LowerBand();
    }

    template <class T2> 
    inline SymBandMatrix(const GenSymBandMatrix<T2>& rhs, int newnlo) : 
      NEW_SIZE(rhs.size(),newnlo)
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      TMVAssert(newnlo <= rhs.nlo());
      if (U==Upper) UpperBand() = rhs.UpperBand().Diags(0,newnlo+1);
      else LowerBand() = rhs.LowerBand().Diags(-newnlo,1);
    }

    template <class T2> 
    inline SymBandMatrix(const GenBandMatrix<T2>& rhs) :
      NEW_SIZE(rhs.rowsize(),U==Upper?rhs.nhi():rhs.nlo())
      { 
#ifdef XTEST_DEBUG
        SetAllTo(T(888));
#endif
        TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
        TMVAssert(rhs.rowsize() == rhs.colsize());
        if (U==Upper) UpperBand() = rhs.UpperBand();
        else LowerBand() = rhs.LowerBand();
      }

    template <class T2> 
    inline SymBandMatrix(const GenMatrix<T2>& rhs, int newnlo) :
      NEW_SIZE(rhs.rowsize(),newnlo)
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      if (U==Upper) UpperBand() = BandMatrixViewOf(rhs,0,nlo());
      else LowerBand() = BandMatrixViewOf(rhs,nlo(),0);
    }

    template <class T2> 
    inline SymBandMatrix(const GenBandMatrix<T2>& rhs, int newnlo) :
      NEW_SIZE(rhs.rowsize(),newnlo)
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      if (U==Upper) UpperBand() = BandMatrixViewOf(rhs,0,nlo());
      else LowerBand() = BandMatrixViewOf(rhs,nlo(),0);
    }

    template <class T2> 
    inline SymBandMatrix(const GenSymMatrix<T2>& rhs, int newnlo) :
      NEW_SIZE(rhs.rowsize(),newnlo)
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      if (U==Upper) UpperBand() = BandMatrixViewOf(rhs.UpperTri(),nlo());
      else LowerBand() = BandMatrixViewOf(rhs.LowerTri(),nlo());
    }

    template <class T2> 
    inline SymBandMatrix(const GenDiagMatrix<T2>& m2) :
      NEW_SIZE(m2.size(),0)
    { 
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      Zero();
      DiagMatrixViewOf(diag()) = m2;
    }

    inline SymBandMatrix(const AssignableToSymBandMatrix<RealType(T)>& m2) :
      NEW_SIZE(m2.size(),m2.nlo())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      TMVAssert(m2.issym());
      m2.AssignTosB(View());
    }

    inline SymBandMatrix(const AssignableToSymBandMatrix<ComplexType(T)>& m2) :
      NEW_SIZE(m2.size(),m2.nlo())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(IsComplex(T()));
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      TMVAssert(m2.issym());
      m2.AssignTosB(View());
    }

    inline SymBandMatrix(const GenDiagMatrix<RealType(T)>& m2) :
      NEW_SIZE(m2.size(),0)
    { 
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      Zero();
      m2.AssignToD(DiagMatrixViewOf(diag()));
      if (nlo() > 0) UpperBandOff().Zero();
    }

    inline SymBandMatrix(const GenDiagMatrix<ComplexType(T)>& m2) :
      NEW_SIZE(m2.size(),0)
    { 
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      TMVAssert(IsComplex(T()));
      Zero();
      m2.AssignToD(DiagMatrixViewOf(diag()));
    }

#undef NEW_SIZE

    virtual inline ~SymBandMatrix() 
    {
#ifdef TMVDEBUG
      SetAllTo(T(999));
#endif
    }

    //
    // Op=
    //

    inline SymBandMatrix<T,U,S,I>& operator=(const SymBandMatrix<T,U,S,I>& m2)
    { 
      TMVAssert(size() == m2.size());
      TMVAssert(nlo() >= m2.nlo());
      if (&m2 != this) {
        if (nlo() > m2.nlo()) LowerBand().Diags(-nlo(),-m2.nlo()).Zero();
        if (S==DiagMajor)
          if (nlo() > m2.nlo())
            std::copy(m2.start_mem(),m2.start_mem()+m2.mem_used(),
                itsm+m2.nlo()*stepi());
          else
            std::copy(m2.start_mem(),m2.start_mem()+linsize,itsm1.get());
        else if (nlo()==m2.nlo())
          std::copy(m2.cptr(),m2.cptr()+m2.mem_used(),itsm1.get());
        else
          LowerBand().Diags(-m2.nlo(),1) = m2.LowerBand();
      }
      return *this;
    }

    template <IndexStyle I2> 
    inline SymBandMatrix<T,U,S,I>& operator=(const SymBandMatrix<T,U,S,I2>& m2)
    { 
      TMVAssert(size() == m2.size());
      TMVAssert(nlo() >= m2.nlo());
      if (nlo() > m2.nlo()) LowerBand().Diags(-nlo(),-m2.nlo()).Zero();
      if (S==DiagMajor)
        if (nlo() > m2.nlo())
          std::copy(m2.start_mem(),m2.start_mem()+m2.mem_used(),
              itsm+m2.nlo()*stepi());
        else
          std::copy(m2.start_mem(),m2.start_mem()+linsize,itsm1.get());
      else if (nlo()==m2.nlo())
        std::copy(m2.cptr(),m2.cptr()+m2.mem_used(),itsm1.get());
      else
        LowerBand().Diags(-m2.nlo(),1) = m2.LowerBand();
      return *this;
    }

    inline SymBandMatrix<T,U,S,I>& operator=(
        const GenSymBandMatrix<RealType(T)>& m2)
    { 
      TMVAssert(size() == m2.size());
      TMVAssert(nlo() >= m2.nlo());
      TMVAssert(m2.issym());
      m2.AssignTosB(View());
      return *this; 
    }

    inline SymBandMatrix<T,U,S,I>& operator=(
        const GenSymBandMatrix<ComplexType(T)>& m2)
    { 
      TMVAssert(IsComplex(T()));
      TMVAssert(size() == m2.size());
      TMVAssert(nlo() >= m2.nlo());
      TMVAssert(m2.issym());
      m2.AssignTosB(View());
      return *this; 
    }

    template <class T2> 
    inline SymBandMatrix<T,U,S,I>& operator=(const GenSymBandMatrix<T2>& m2)
    { 
      TMVAssert(size() == m2.size());
      TMVAssert(IsReal(T2()) || m2.issym());
      LowerBand() = m2.LowerBand();
      return *this; 
    }

    inline SymBandMatrix<T,U,S,I>& operator=(T x) 
    { return SetToIdentity(x); }

    inline SymBandMatrix<T,U,S,I>& operator=(
        const AssignableToSymBandMatrix<RealType(T)>& m2)
    { 
      TMVAssert(size() == m2.size());
      TMVAssert(m2.issym());
      m2.AssignTosB(View());
      return *this;
    }

    inline SymBandMatrix<T,U,S,I>& operator=(
        const AssignableToSymBandMatrix<ComplexType(T)>& m2)
    { 
      TMVAssert(IsComplex(T()));
      TMVAssert(size() == m2.size());
      TMVAssert(m2.issym());
      m2.AssignTosB(View());
      return *this;
    }

    inline SymBandMatrix<T,U,S,I>& operator=(
        const GenDiagMatrix<RealType(T)>& m2)
    { 
      TMVAssert(size() == m2.size());
      View() = m2;
      return *this;
    }

    inline SymBandMatrix<T,U,S,I>& operator=(
        const GenDiagMatrix<ComplexType(T)>& m2) 
    { 
      TMVAssert(IsComplex(T()));
      TMVAssert(size() == m2.size());
      TMVAssert(issym());
      View() = m2;
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
        return okij(i,j) ? cref(i,j) : T(0);
      } else {
        TMVAssert(i>0 && i<=int(size()));
        TMVAssert(j>0 && j<=int(size()));
        return okij(i-1,j-1) ? cref(i-1,j-1) : T(0);
      }
    }

    inline T& operator()(int i, int j) 
    { 
      if (I==CStyle) {
        TMVAssert(i>=0 && i<int(size()));
        TMVAssert(j>=0 && j<int(size()));
        TMVAssert(okij(i,j));
        return ref(i,j); 
      } else {
        TMVAssert(i>0 && i<=int(size()));
        TMVAssert(j>0 && j<=int(size()));
        TMVAssert(okij(i-1,j-1));
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
      TMVAssert(j1==j2 || okij(i,j1));
      TMVAssert(j1==j2 || okij(i,j2-1));
      TMVAssert( i-j1<=0 || j2-i<=1 );
      if ((U==Upper && i-j1<=0) || (U==Lower && j2-i<=1))
        return ConstVectorView<T,I>(itsm+i*stepi()+j1*stepj(),
            j2-j1,stepj(),NonConj); 
      else
        return ConstVectorView<T,I>(itsm+i*stepj()+j1*stepi(),
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
      TMVAssert(i1==i2 || okij(i1,j));
      TMVAssert(i1==i2 || okij(i2-1,j));
      TMVAssert( j-i1<=0 || i2-j<=1 );
      if ((U==Upper && i2-j<=1) || (U==Lower && j-i1<=0))
        return ConstVectorView<T,I>(itsm+i1*stepi()+j*stepj(),
            i2-i1,stepi(),NonConj); 
      else 
        return ConstVectorView<T,I>(itsm+i1*stepj()+j*stepi(),
            i2-i1,stepj(),NonConj); 
    }

    inline ConstVectorView<T,I> diag() const
    { return ConstVectorView<T,I>(itsm,size(),diagstep(),NonConj); }

    inline ConstVectorView<T,I> diag(int i) const
    {
      TMVAssert(i>=-nlo() && i<=nlo());
      TMVAssert(i>=-int(size()) && i<=int(size())); 
      i = std::abs(i);
      return ConstVectorView<T,I>(itsm+i*(U==Upper?stepj():stepi()),
          size()-i,diagstep(),NonConj); 
    }

    inline ConstVectorView<T,I> diag(int i, int j1, int j2) const
    {
      TMVAssert(i>=-nlo() && i<=nlo());
      TMVAssert(i>=-int(size()) && i<=int(size())); 
      if (I==FortranStyle) { 
        TMVAssert(j1>0 && j1-j2<=0 && j2<=int(size())-std::abs(i)); --j1;
      } else {
        TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size())-std::abs(i));
      }
      i = std::abs(i);
      return ConstVectorView<T,I>(
          itsm+i*(U==Upper?stepj():stepi())+j1*diagstep(),
          j2-j1,diagstep(),NonConj);
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
      TMVAssert(j1==j2 || okij(i,j1));
      TMVAssert(j1==j2 || okij(i,j2-1));
      TMVAssert( i-j1<=0 || j2-i<=1 );
      if ((U==Upper && i-j1<=0) || (U==Lower && j2-i<=1))
        return VectorView<T,I>(itsm+i*stepi()+j1*stepj(),
            j2-j1,stepj(),NonConj FIRSTLAST); 
      else
        return VectorView<T,I>(itsm+i*stepj()+j1*stepi(),
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
      TMVAssert(i1==i2 || okij(i1,j));
      TMVAssert(i1==i2 || okij(i2-1,j));
      TMVAssert( j-i1<=0 || i2-j<=1 );
      if ((U==Upper && i2-j<=1) || (U==Lower && j-i1<=0))
        return VectorView<T,I>(itsm+i1*stepi()+j*stepj(),
            i2-i1,stepi(),NonConj FIRSTLAST); 
      else 
        return VectorView<T,I>(itsm+i1*stepj()+j*stepi(),
            i2-i1,stepj(),NonConj FIRSTLAST); 
    }

    inline VectorView<T,I> diag()
    { return VectorView<T,I>(itsm,size(),diagstep(),NonConj FIRSTLAST); }

    inline VectorView<T,I> diag(int i) 
    {
      TMVAssert(i>=-nlo() && i<=nlo());
      TMVAssert(i>=-int(size()) && i<=int(size())); 
      i = std::abs(i);
      return VectorView<T,I>(itsm+i*(U==Upper?stepj():stepi()),
          size()-i,diagstep(),NonConj FIRSTLAST); 
    }

    inline VectorView<T,I> diag(int i, int j1, int j2) 
    {
      TMVAssert(i>=-nlo() && i<=nlo());
      TMVAssert(i>=-int(size()) && i<=int(size())); 
      if (I==FortranStyle) { 
        TMVAssert(j1>0 && j1-j2<=0 && j2<=int(size())-std::abs(i)); --j1;
      } else {
        TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size())-std::abs(i));
      }
      i = std::abs(i);
      return VectorView<T,I>(
          itsm+i*(U==Upper?stepj():stepi())+j1*diagstep(),
          j2-j1,diagstep(),NonConj FIRSTLAST);
    }


    //
    // Modifying Functions
    //

    inline SymBandMatrix<T,U,S,I>& Zero() 
    {
      std::fill_n(itsm1.get(),linsize,T(0));
      return *this;
    }

    inline SymBandMatrix<T,U,S,I>& SetAllTo(T x) 
    { UpperBand().SetAllTo(x); return *this; }

    inline SymBandMatrix<T,U,S,I>& Clip(RealType(T) thresh) 
    { UpperBand().Clip(thresh); return *this; }

    inline SymBandMatrix<T,U,S,I>& ConjugateSelf() 
    { if (IsComplex(T())) UpperBand().ConjugateSelf(); return *this; }

    inline SymBandMatrix<T,U,S,I>& TransposeSelf() 
    { return *this; }

    inline SymBandMatrix<T,U,S,I>& SetToIdentity(T x=T(1)) 
    { Zero(); diag().SetAllTo(x); return *this; }


    //
    // SubMatrix
    //

    inline ConstMatrixView<T,I> SubMatrix(int i1, int i2, int j1, int j2) const
    {
      TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,1,1));
      if (I==FortranStyle) { --i1; --j1; }
      TMVAssert(i2-1<=j1 || j2-1<=i1);
      if ((U==Upper && i2-1<=j1) || (U==Lower && j2-1<=i1))
        return ConstMatrixView<T,I>(itsm+i1*stepi()+j1*stepj(),
            i2-i1,j2-j1,stepi(),stepj(),S,NonConj);
      else
        return ConstMatrixView<T,I>(itsm+i1*stepj()+j1*stepi(),
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
        return ConstMatrixView<T,I>(itsm+i1*stepi()+j1*stepj(),
            (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
            newstor,NonConj);
      else
        return ConstMatrixView<T,I>(itsm+i1*stepj()+j1*stepi(),
            (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
            TransOf(newstor),NonConj);
    }

    inline ConstBandMatrixView<T,I> SubBandMatrix(
        int i1, int i2, int j1, int j2, int newnlo, int newnhi) const
    {
      TMVAssert(View().OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
      if (I==FortranStyle) { --i1; --j1; }
      TMVAssert(i1+newnlo-j1<=0 || j1+newnhi-i1<=0);
      if ((i1+newnlo-j1<=0 && U==Upper) || (j1+newnhi-i1<=0 && U==Lower))
        return ConstBandMatrixView<T,I>(itsm+i1*stepi()+j1*stepj(),
            i2-i1,j2-j1,newnlo,newnhi,stepi(),stepj(),diagstep(),
            S,NonConj);
      else
        return ConstBandMatrixView<T,I>(itsm+i1*stepj()+j1*stepi(),
            i2-i1,j2-j1,newnlo,newnhi,stepj(),stepi(),diagstep(),
            TransOf(S),NonConj);
    }

    inline ConstBandMatrixView<T,I> SubBandMatrix(
        int i1, int i2, int j1, int j2, int newnlo, int newnhi,
        int istep, int jstep) const
    {
      TMVAssert(View().OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,
            istep,jstep));
      if (I==FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
      const StorageType newstor =
      iscm() ? (istep == 1 ? ColMajor : NoMajor) :
      isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
      TMVAssert(i1+newnlo*istep<=j1 || j1+newnhi*jstep<=i1);
      if ((i1+newnlo*istep<=j1 && U==Upper) || 
          (j1+newnhi*jstep<=i1 && U==Lower)) {
        const int newstepi = stepi()*istep;
        const int newstepj = stepj()*jstep;
        return ConstBandMatrixView<T,I>(itsm+i1*stepi()+j1*stepj(),
            (i2-i1)/istep,(j2-j1)/jstep,newnlo,newnhi,
            newstepi,newstepj,newstepi+newstepj,newstor,NonConj);
      } else {
        const int newstepi = stepj()*istep;
        const int newstepj = stepi()*jstep;
        return ConstBandMatrixView<T,I>(itsm+i1*stepj()+j1*stepi(),
            (i2-i1)/istep,(j2-j1)/jstep,newnlo,newnhi,
            newstepi,newstepj,newstepi+newstepj,
            TransOf(newstor),NonConj);
      }
    }

    inline ConstVectorView<T,I> SubVector(int i, int j,
        int istep, int jstep, int n) const
    {
      TMVAssert(View().OKSubVector(i,j,istep,jstep,n));
      if (I==FortranStyle) { --i; --j; }
      if ((U==Upper && i-j<=0) || (U==Lower && j-i<=0))
        return ConstVectorView<T,I>(itsm+i*stepi()+j*stepj(),n,
            istep*stepi()+jstep*stepj(),NonConj);
      else
        return ConstVectorView<T,I>(itsm+i*stepj()+j*stepi(),n,
            istep*stepj()+jstep*stepi(),NonConj);
    }

    inline ConstSymMatrixView<T,I> SubSymMatrix(int i1, int i2) const
    {
      TMVAssert(View().OKSubSymMatrix(i1,i2,1));
      if (I==FortranStyle) { --i1; }
      return ConstSymMatrixView<T,I>(itsm+i1*diagstep(),i2-i1,
          stepi(),stepj(),Sym,U,S,NonConj);
    }

    inline ConstSymMatrixView<T,I> SubSymMatrix(int i1, int i2,
        int istep) const
    {
      TMVAssert(View().OKSubSymMatrix(i1,i2,istep));
      if (I==FortranStyle) { --i1; i2+=istep-1; }
      return ConstSymMatrixView<T,I>(itsm+i1*diagstep(),
          (i2-i1)/istep,istep*stepi(),istep*stepj(),Sym,U,
          istep==1 ? S : NoMajor, NonConj);
    }

    inline ConstSymBandMatrixView<T,I> SubSymBandMatrix(
        int i1, int i2, int newnlo) const
    {
      TMVAssert(View().OKSubSymBandMatrix(i1,i2,newnlo,1));
      if (I==FortranStyle) { --i1; }
      return ConstSymBandMatrixView<T,I>(itsm+i1*diagstep(),i2-i1,
          newnlo,stepi(),stepj(),diagstep(),Sym,U,S,NonConj);
    }

    inline ConstSymBandMatrixView<T,I> SubSymBandMatrix(int i1, int i2) const
    { 
      return SubSymBandMatrix(i1,i2,MIN(nlo(),i2-i1-(I==CStyle?0:1))); 
    }

    inline ConstSymBandMatrixView<T,I> SubSymBandMatrix(
        int i1, int i2, int newnlo, int istep) const
    {
      TMVAssert(View().OKSubSymBandMatrix(i1,i2,newnlo,istep));
      if (I==FortranStyle) { --i1; i2+=istep-1; }
      return ConstSymBandMatrixView<T,I>(
          itsm+i1*diagstep(), (i2-i1)/istep, newnlo,
          istep*stepi(),istep*stepj(),istep*diagstep(),
          Sym,U,istep==1 ? S : NoMajor,NonConj);
    }

    inline ConstSymBandMatrixView<T,I> SymDiags(int newnlo) const
    {
      TMVAssert(newnlo>=0 && newnlo <= nlo());
      return ConstSymBandMatrixView<T,I>(itsm,size(),newnlo,
          stepi(),stepj(),diagstep(),Sym,U,S,NonConj);
    }

    inline ConstBandMatrixView<T> Diags(int k1, int k2) const
    {
      TMVAssert(k1>=-nlo() && k1<k2 && k2<=nlo()+1);
      TMVAssert(k2 <= 1 || k1 >= 0);

      if (k1 >= 0) {
        const int newsize = int(size())-k1;
        const int newnhi = k2-k1-1;
        if (uplo() == Upper)
          return ConstBandMatrixView<T>(itsm+k1*stepj(),
              newsize, newsize, 0, newnhi, stepi(), stepj(), diagstep(),
              stor(), ct());
        else
          return ConstBandMatrixView<T>(itsm+k1*stepi(),
              newsize, newsize, 0, newnhi, stepj(), stepi(), diagstep(),
              TransOf(stor()), issym()?ct():ConjOf(T,ct()));
      } else {
        const int newsize = int(size())+k2-1;
        const int newnlo = k2-k1-1;
        if (uplo() == Lower)
          return ConstBandMatrixView<T>(itsm-k2*stepi(),
              newsize, newsize, newnlo, 0, stepi(), stepj(), diagstep(),
              stor(), ct());
        else
          return ConstBandMatrixView<T>(itsm-k2*stepj(),
              newsize, newsize, newnlo, 0, stepj(), stepi(), diagstep(),
              TransOf(stor()), issym()?ct():ConjOf(T,ct()));
      }
    }

    inline ConstBandMatrixView<T> UpperBand() const
    {
      if (uplo() == Upper)
        return ConstBandMatrixView<T>(itsm,size(),size(),0,nlo(),
            stepi(),stepj(),diagstep(),S,NonConj);
      else
        return ConstBandMatrixView<T>(itsm,size(),size(),0,nlo(),
            stepj(),stepi(),diagstep(),TransOf(S),NonConj);
    }

    inline ConstBandMatrixView<T> UpperBandOff() const
    {
      if (uplo() == Upper)
        return ConstBandMatrixView<T>(itsm+stepj(),size()-1,size()-1,
            0,nlo()-1,stepi(),stepj(),diagstep(),S,NonConj);
      else
        return ConstBandMatrixView<T>(itsm+stepi(),size()-1,size()-1,
            0,nlo()-1,stepj(),stepi(),diagstep(),TransOf(S),NonConj);
    }

    inline ConstBandMatrixView<T> LowerBand() const
    {
      if (uplo() == Lower)
        return ConstBandMatrixView<T>(itsm,size(),size(),nlo(),0,
            stepi(),stepj(),diagstep(),S,NonConj);
      else
        return ConstBandMatrixView<T>(itsm,size(),size(),nlo(),0,
            stepj(),stepi(),diagstep(),TransOf(S),NonConj);
    }

    inline ConstBandMatrixView<T> LowerBandOff() const
    {
      if (uplo() == Lower)
        return ConstBandMatrixView<T>(itsm+stepi(),size()-1,size()-1,
            nlo()-1,0,stepi(),stepj(),diagstep(),S,NonConj);
      else
        return ConstBandMatrixView<T>(itsm+stepj(),size()-1,size()-1,
            nlo()-1,0,stepj(),stepi(),diagstep(),TransOf(S),NonConj);
    }

    inline MatrixView<T,I> SubMatrix(int i1, int i2, int j1, int j2)
    {
      TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,1,1));
      if (I==FortranStyle) { --i1; --j1; }
      TMVAssert(i2-1<=j1 || j2-1<=i1);
      if ((U==Upper && i2-1<=j1) || (U==Lower && j2-1<=i1))
        return MatrixView<T,I>(itsm+i1*stepi()+j1*stepj(),
            i2-i1,j2-j1,stepi(),stepj(),S,NonConj FIRSTLAST);
      else
        return MatrixView<T,I>(itsm+i1*stepj()+j1*stepi(),
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
        return MatrixView<T,I>(itsm+i1*stepi()+j1*stepj(),
            (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
            newstor,NonConj FIRSTLAST);
      else
        return MatrixView<T,I>(itsm+i1*stepj()+j1*stepi(),
            (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
            TransOf(newstor),NonConj FIRSTLAST);
    }

    inline BandMatrixView<T,I> SubBandMatrix(
        int i1, int i2, int j1, int j2, int newnlo, int newnhi)
    {
      TMVAssert(View().OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
      if (I==FortranStyle) { --i1; --j1; }
      TMVAssert(i1+newnlo-j1<=0 || j1+newnhi-i1<=0);
      if ((i1+newnlo-j1<=0 && U==Upper) || (j1+newnhi-i1<=0 && U==Lower))
        return BandMatrixView<T,I>(itsm+i1*stepi()+j1*stepj(),
            i2-i1,j2-j1,newnlo,newnhi,stepi(),stepj(),diagstep(),
            S,NonConj FIRSTLAST);
      else
        return BandMatrixView<T,I>(itsm+i1*stepj()+j1*stepi(),
            i2-i1,j2-j1,newnlo,newnhi,stepj(),stepi(),diagstep(),
            TransOf(S),NonConj FIRSTLAST);
    }

    inline BandMatrixView<T,I> SubBandMatrix(
        int i1, int i2, int j1, int j2, int newnlo, int newnhi,
        int istep, int jstep)
    {
      TMVAssert(View().OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,
            istep,jstep));
      if (I==FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
      const StorageType newstor=
      iscm() ? (istep == 1 ? ColMajor : NoMajor) :
      isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
      TMVAssert(i1+newnlo*istep<=j1 || j1+newnhi*jstep<=i1);
      if ((i1+newnlo*istep<=j1 && U==Upper) || 
          (j1+newnhi*jstep<=i1 && U==Lower)) {
        const int newstepi = stepi()*istep;
        const int newstepj = stepj()*jstep;
        return BandMatrixView<T,I>(itsm+i1*stepi()+j1*stepj(),
            (i2-i1)/istep,(j2-j1)/jstep,newnlo,newnhi,
            newstepi,newstepj,newstepi+newstepj,newstor,NonConj FIRSTLAST);
      } else {
        const int newstepi = stepj()*istep;
        const int newstepj = stepi()*jstep;
        return BandMatrixView<T,I>(itsm+i1*stepj()+j1*stepi(),
            (i2-i1)/istep,(j2-j1)/jstep,newnlo,newnhi,
            newstepi,newstepj,newstepi+newstepj,
            TransOf(newstor),NonConj FIRSTLAST);
      }
    }

    inline VectorView<T,I> SubVector(int i, int j,
        int istep, int jstep, int n)
    {
      TMVAssert(View().OKSubVector(i,j,istep,jstep,n));
      if (I==FortranStyle) { --i; --j; }
      if ((U==Upper && i-j<=0) || (U==Lower && j-i<=0))
        return VectorView<T,I>(itsm+i*stepi()+j*stepj(),n,
            istep*stepi()+jstep*stepj(),NonConj FIRSTLAST);
      else
        return VectorView<T,I>(itsm+i*stepj()+j*stepi(),n,
            istep*stepj()+jstep*stepi(),NonConj FIRSTLAST);
    }

    inline SymMatrixView<T,I> SubSymMatrix(int i1, int i2)
    {
      TMVAssert(View().OKSubSymMatrix(i1,i2,1));
      if (I==FortranStyle) { --i1; }
      return SymMatrixView<T,I>(itsm+i1*diagstep(),i2-i1,
          stepi(),stepj(),Sym,U,S,NonConj FIRSTLAST);
    }

    inline SymMatrixView<T,I> SubSymMatrix(int i1, int i2, int istep) 
    {
      TMVAssert(View().OKSubSymMatrix(i1,i2,istep));
      if (I==FortranStyle) { --i1; i2+=istep-1; }
      return SymMatrixView<T,I>(itsm+i1*diagstep(),(i2-i1)/istep,
          istep*stepi(),istep*stepj(),Sym,U,
          istep==1 ? S : NoMajor,NonConj FIRSTLAST);
    }

    inline SymBandMatrixView<T,I> SubSymBandMatrix(int i1, int i2, int newnlo)
    {
      TMVAssert(View().OKSubSymBandMatrix(i1,i2,newnlo,1));
      if (I==FortranStyle) { --i1; }
      return SymBandMatrixView<T,I>(itsm+i1*diagstep(),i2-i1,newnlo,
          stepi(),stepj(),diagstep(),Sym,U,S,NonConj FIRSTLAST);
    }

    inline SymBandMatrixView<T,I> SubSymBandMatrix(int i1, int i2)
    { 
      return SubSymBandMatrix(i1,i2,MIN(nlo(),i2-i1-(I==CStyle?0:1))); 
    }

    inline SymBandMatrixView<T,I> SubSymBandMatrix(
        int i1, int i2, int newnlo, int istep)
    {
      TMVAssert(View().OKSubSymBandMatrix(i1,i2,newnlo,istep));
      if (I==FortranStyle) { --i1; i2+=istep-1; }
      return SymBandMatrixView<T,I>(itsm+i1*diagstep(),(i2-i1)/istep,
          newnlo,istep*stepi(),istep*stepj(),istep*diagstep(),
          Sym,U,istep==1 ? S : NoMajor,NonConj FIRSTLAST);
    }

    inline SymBandMatrixView<T,I> SymDiags(int newnlo)
    {
      TMVAssert(newnlo>=0 && newnlo <= nlo());
      return SymBandMatrixView<T,I>(itsm,size(),newnlo,
          stepi(),stepj(),diagstep(),Sym,U,S,NonConj FIRSTLAST);
    }

    inline BandMatrixView<T> Diags(int k1, int k2) 
    {
      TMVAssert(k1>=-nlo() && k1<k2 && k2<=nlo()+1);
      TMVAssert(k2 <= 1 || k1 >= 0);

      if (k1 >= 0) {
        const int newsize = int(size())-k1;
        const int newnhi = k2-k1-1;
        if (uplo() == Upper)
          return BandMatrixView<T>(itsm+k1*stepj(),
              newsize, newsize, 0, newnhi, stepi(), stepj(), diagstep(),
              stor(), ct() FIRSTLAST );
        else
          return BandMatrixView<T>(itsm+k1*stepi(),
              newsize, newsize, 0, newnhi, stepj(), stepi(), diagstep(),
              TransOf(stor()), issym()?ct():ConjOf(T,ct()) FIRSTLAST );
      } else {
        const int newsize = int(size())+k2-1;
        const int newnlo = k2-k1-1;
        if (uplo() == Lower)
          return BandMatrixView<T>(itsm-k2*stepi(),
              newsize, newsize, newnlo, 0, stepi(), stepj(), diagstep(),
              stor(), ct() FIRSTLAST );
        else
          return BandMatrixView<T>(itsm-k2*stepj(),
              newsize, newsize, newnlo, 0, stepj(), stepi(), diagstep(),
              TransOf(stor()), issym()?ct():ConjOf(T,ct()) FIRSTLAST );
      }
    }

    inline BandMatrixView<T> UpperBand()
    {
      if (uplo() == Upper)
        return BandMatrixView<T>(itsm,size(),size(),0,nlo(),
            stepi(),stepj(),diagstep(),S,NonConj FIRSTLAST);
      else
        return BandMatrixView<T>(itsm,size(),size(),0,nlo(),
            stepj(),stepi(),diagstep(),TransOf(S),NonConj FIRSTLAST);
    }

    inline BandMatrixView<T> UpperBandOff()
    {
      if (uplo() == Upper)
        return BandMatrixView<T>(itsm+stepj(),size()-1,size()-1,
            0,nlo()-1,stepi(),stepj(),diagstep(),S,NonConj FIRSTLAST);
      else
        return BandMatrixView<T>(itsm+stepi(),size()-1,size()-1,
            0,nlo()-1,stepj(),stepi(),diagstep(),TransOf(S),NonConj 
            FIRSTLAST);
    }

    inline BandMatrixView<T> LowerBand()
    {
      if (uplo() == Lower)
        return BandMatrixView<T>(itsm,size(),size(),nlo(),0,
            stepi(),stepj(),diagstep(),S,NonConj FIRSTLAST);
      else
        return BandMatrixView<T>(itsm,size(),size(),nlo(),0,
            stepj(),stepi(),diagstep(),TransOf(S),NonConj FIRSTLAST);
    }

    inline BandMatrixView<T> LowerBandOff()
    {
      if (uplo() == Lower)
        return BandMatrixView<T>(itsm+stepi(),size()-1,size()-1,
            nlo()-1,0,stepi(),stepj(),diagstep(),S,NonConj FIRSTLAST);
      else
        return BandMatrixView<T>(itsm+stepj(),size()-1,size()-1,
            nlo()-1,0,stepj(),stepi(),diagstep(),TransOf(S),NonConj 
            FIRSTLAST);
    }

    inline ConstSymBandMatrixView<T,I> View() const
    { 
      return ConstSymBandMatrixView<T,I>(itsm,size(),nlo(),
          stepi(),stepj(),diagstep(),Sym,U,S,NonConj);
    }

    inline ConstSymBandMatrixView<T,I> Transpose() const
    {
      return ConstSymBandMatrixView<T,I>(itsm,size(),nlo(),
          stepj(),stepi(),diagstep(),Sym,UTransOf(U),TransOf(S),NonConj);
    }

    inline ConstSymBandMatrixView<T,I> Conjugate() const
    { 
      return ConstSymBandMatrixView<T,I>(itsm,size(),nlo(),
          stepi(),stepj(),diagstep(),Sym,U,S,ConjOf(T,NonConj));
    }

    inline ConstSymBandMatrixView<T,I> Adjoint() const
    {
      return ConstSymBandMatrixView<T,I>(itsm,size(),nlo(),
          stepj(),stepi(),diagstep(),
          Sym,UTransOf(U),TransOf(S),ConjOf(T,NonConj));
    }

    inline ConstSymBandMatrixView<RealType(T),I> Real() const
    {
      return ConstSymBandMatrixView<RealType(T)>(
          reinterpret_cast<const RealType(T)*>(itsm),size(),nlo(),
          IsReal(T()) ? stepi() : 2*stepi(),
          IsReal(T()) ? stepj() : 2*stepj(),
          IsReal(T()) ? diagstep() : 2*diagstep(),
          Sym,U,IsReal(T())?S:NoMajor,NonConj);
    }

    inline ConstSymBandMatrixView<RealType(T),I> Imag() const
    {
      TMVAssert(IsComplex(T()));
      return ConstSymBandMatrixView<RealType(T)>(
          reinterpret_cast<const RealType(T)*>(itsm)+1,size(),nlo(),
          2*stepi(),2*stepj(),2*diagstep(),Sym,U,NoMajor,NonConj);
    } 

    inline SymBandMatrixView<T,I> View() 
    { 
      return SymBandMatrixView<T,I>(itsm,size(),nlo(),
          stepi(),stepj(),diagstep(), Sym,U,S,NonConj FIRSTLAST);
    }

    inline SymBandMatrixView<T,I> Transpose() 
    {
      return SymBandMatrixView<T,I>(itsm,size(),nlo(),
          stepj(),stepi(),diagstep(),
          Sym,UTransOf(U),TransOf(S),NonConj FIRSTLAST);
    }

    inline SymBandMatrixView<T,I> Conjugate() 
    { 
      return SymBandMatrixView<T,I>(itsm,size(),nlo(),
          stepi(),stepj(),diagstep(),
          Sym,U,S,ConjOf(T,NonConj) FIRSTLAST);
    }

    inline SymBandMatrixView<T,I> Adjoint() 
    {
      return SymBandMatrixView<T,I>(itsm,size(),nlo(),
          stepj(),stepi(),diagstep(),
          Sym,UTransOf(U),TransOf(S),ConjOf(T,NonConj) FIRSTLAST);
    }

    inline SymBandMatrixView<RealType(T),I> Real()
    {
      return SymBandMatrixView<RealType(T)>(
          reinterpret_cast<RealType(T)*>(itsm),size(),nlo(),
          IsReal(T()) ? stepi() : 2*stepi(),
          IsReal(T()) ? stepj() : 2*stepj(),
          IsReal(T()) ? diagstep() : 2*diagstep(),
          Sym,U,IsReal(T())?S:NoMajor,NonConj
#ifdef TMVFLDEBUG
          ,reinterpret_cast<const RealType(T)*>(first)
          ,reinterpret_cast<const RealType(T)*>(last)
#endif
          );
    }

    inline SymBandMatrixView<RealType(T),I> Imag()
    {
      TMVAssert(IsComplex(T()));
      return SymBandMatrixView<RealType(T)>(
          reinterpret_cast<RealType(T)*>(itsm)+1,size(),nlo(),
          2*stepi(),2*stepj(),2*diagstep(),Sym,U,NoMajor,NonConj
#ifdef TMVFLDEBUG
          ,reinterpret_cast<const RealType(T)*>(first)+1
          ,reinterpret_cast<const RealType(T)*>(last)+1
#endif
          );
    } 

    inline size_t size() const { return itss; }
    inline int nlo() const { return itslo; }
    inline size_t mem_used() const { return linsize; }
    inline const T* start_mem() const { return itsm1.get(); }
    inline const T* cptr() const { return itsm; }
    inline T* ptr() { return itsm; }
    inline int stepi() const { return itssi; }
    inline int stepj() const { return itssj; }
    inline int diagstep() const { return itssd; }
    inline SymType sym() const { return Sym; }
    inline UpLoType uplo() const { return U; }
    inline StorageType stor() const { return S; }
    inline ConjType ct() const { return NonConj; }
    inline bool isrm() const { return S==RowMajor; }
    inline bool iscm() const { return S==ColMajor; }
    inline bool isdm() const { return S==DiagMajor; }
    inline bool isconj() const { return false; }
    inline bool isherm() const { return IsReal(T()); }
    inline bool issym() const { return true; }

    inline T& ref(int i, int j)
    {
      if ((U==Upper && i <= j) || (U==Lower && i>=j)) 
        return itsm[i*itssi + j*itssj];
      else 
        return itsm[j*itssi + i*itssj];
    }

    inline T cref(int i, int j) const 
    {
      if ((U==Upper && i <= j) || (U==Lower && i>=j)) 
        return itsm[i*itssi + j*itssj];
      else 
        return itsm[j*itssi + i*itssj];
    }

  protected :

    const size_t linsize;
    auto_array<T> itsm1;
    const size_t itss;
    const int itslo;
    const int itssi;
    const int itssj;
    const int itssd;
    T*const itsm;

#ifdef TMVFLDEBUG
  public:
    const T*const first;
    const T*const last;
  protected :
#endif

    inline bool okij(int i, int j) const
    { return (j+nlo() >= i && i+nlo() >= j); }

  }; // SymBandMatrix

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  class HermBandMatrix : 
    public GenSymBandMatrix<T>
  {

  public:

    //
    // Constructors
    //

#define NEW_SIZE(s,lo) \
    linsize(BandStorageLength(S,s,s,lo,0)), \
    itsm1(new T[linsize]), itss(s), itslo(lo), \
    itssi(S==DiagMajor ? -int(s)+1 : S==RowMajor ? lo : 1), \
    itssj(S==DiagMajor ? int(s) : S==RowMajor ? 1 : lo), \
    itssd(S==DiagMajor ? 1 : lo+1), \
    itsm((S==DiagMajor && U==Lower) ? itsm1.get()-lo*itssi : itsm1.get()) \
    DEFFIRSTLAST(itsm1.get(),itsm1.get()+linsize) 

    inline HermBandMatrix(size_t s, int lo) : NEW_SIZE(s,lo) 
    { 
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
#ifdef TMVDEBUG
      SetAllTo(T(888));
#endif
    }

    inline HermBandMatrix(size_t s, int lo, const RealType(T) x) : 
      NEW_SIZE(s,lo) 
    { 
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      SetAllTo(x);
    }

    inline HermBandMatrix(size_t s, int lo, const T* vv) : NEW_SIZE(s,lo)
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      std::copy(vv,vv+linsize,itsm1.get());
      TMVAssert(this->HermOK());
    }

    inline HermBandMatrix(size_t s, int lo, const std::vector<T>& vv) : 
      NEW_SIZE(s,lo)
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      TMVAssert(vv.size() == linsize);
      T* vi = itsm1.get();
      typename std::vector<T>::const_iterator vvi = vv.begin();
      for(int i=linsize;i>0;--i,++vi,++vvi) *vi = *vvi;
      TMVAssert(this->HermOK());
    }

    inline HermBandMatrix(const HermBandMatrix<T,U,S,I>& rhs) : 
      NEW_SIZE(rhs.size(),rhs.nlo())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      std::copy(rhs.start_mem(),rhs.start_mem()+linsize,itsm1.get());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    template <UpLoType U2, StorageType S2, IndexStyle I2> 
    inline HermBandMatrix(const HermBandMatrix<T,U2,S2,I2>& rhs) : 
      NEW_SIZE(rhs.size(),rhs.nlo())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      if (U==U2 && S==S2)
        std::copy(rhs.start_mem(),rhs.start_mem()+linsize,itsm1.get());
      else if (S!=DiagMajor && (S2==RowMajor || S2==ColMajor) && 
          S!=S2 && U!=U2) {
        std::copy(rhs.start_mem(),rhs.start_mem()+linsize,itsm1.get());
        ConjugateSelf();
      }
      else if (U==Upper) UpperBand() = rhs.UpperBand();
      else LowerBand() = rhs.LowerBand();
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    inline HermBandMatrix(const GenSymBandMatrix<RealType(T)>& rhs) :
      NEW_SIZE(rhs.size(),rhs.nlo())
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      if (rhs.isherm()) {
        rhs.AssignTosB(View());
      } else {
        if (U==Upper) 
          UpperBand() = rhs.UpperBand();
        else 
          LowerBand() = rhs.LowerBand();
      }
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    inline HermBandMatrix(const GenSymBandMatrix<ComplexType(T)>& rhs) :
      NEW_SIZE(rhs.size(),rhs.nlo())
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(IsComplex(T()));
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      if (rhs.isherm()) {
        rhs.AssignTosB(View());
      } else {
        if (U==Upper) 
          UpperBand() = rhs.UpperBand();
        else 
          LowerBand() = rhs.LowerBand();
        if (IsComplex(T())) diag().Imag().Zero();
      }
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    template <class T2> 
    inline HermBandMatrix(const GenSymBandMatrix<T2>& rhs) :
      NEW_SIZE(rhs.size(),rhs.nlo())
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      if (U==Upper) 
        UpperBand() = rhs.UpperBand();
      else 
        LowerBand() = rhs.LowerBand();
      if (IsComplex(T())) diag().Imag().Zero();
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    template <class T2> 
    inline HermBandMatrix(const GenSymBandMatrix<T2>& rhs, int newnlo) : 
      NEW_SIZE(rhs.size(),newnlo)
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      TMVAssert(newnlo <= rhs.nlo());
      if (U==Upper) UpperBand() = rhs.UpperBand().Diags(0,newnlo+1);
      else LowerBand() = rhs.LowerBand().Diags(-newnlo,1);
      if (IsComplex(T())) diag().Imag().Zero();
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    template <class T2> 
    inline HermBandMatrix(const GenBandMatrix<T2>& rhs) :
      NEW_SIZE(rhs.rowsize(),rhs.nlo())
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      if (U==Upper) UpperBand() = rhs.UpperBand();
      else LowerBand() = rhs.LowerBand();
      if (IsComplex(T())) diag().Imag().Zero();
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    template <class T2> 
    inline HermBandMatrix(const GenMatrix<T2>& rhs,
        int newnlo) :
      NEW_SIZE(rhs.rowsize(),newnlo)
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      if (U==Upper) UpperBand() = BandMatrixViewOf(rhs,0,nlo());
      else LowerBand() = BandMatrixViewOf(rhs,nlo(),0);
      if (IsComplex(T())) diag().Imag().Zero();
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    template <class T2> 
    inline HermBandMatrix(const GenBandMatrix<T2>& rhs, int newnlo) :
      NEW_SIZE(rhs.rowsize(),newnlo)
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      if (U==Upper) UpperBand() = BandMatrixViewOf(rhs,0,nlo());
      else LowerBand() = BandMatrixViewOf(rhs,nlo(),0);
      if (IsComplex(T())) diag().Imag().Zero();
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    template <class T2> 
    inline HermBandMatrix(const GenSymMatrix<T2>& rhs,
        int newnlo) :
      NEW_SIZE(rhs.rowsize(),newnlo)
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      if (U==Upper) 
        UpperBand() = BandMatrixViewOf(rhs.UpperTri(),nlo());
      else 
        LowerBand() = BandMatrixViewOf(rhs.LowerTri(),nlo());
      if (IsComplex(T())) diag().Imag().Zero();
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    template <class T2> 
    inline HermBandMatrix(const GenDiagMatrix<T2>& m2) :
      NEW_SIZE(m2.size(),0)
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      Zero();
      DiagMatrixViewOf(diag()) = m2;
      if (IsComplex(T())) diag().Imag().Zero();
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    inline HermBandMatrix(const AssignableToSymBandMatrix<RealType(T)>& m2) :
      NEW_SIZE(m2.size(),m2.nlo())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      TMVAssert(m2.isherm());
      m2.AssignTosB(View());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    inline HermBandMatrix(const AssignableToSymBandMatrix<ComplexType(T)>& m2) :
      NEW_SIZE(m2.size(),m2.nlo())
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(IsComplex(T()));
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      TMVAssert(m2.isherm());
      m2.AssignTosB(View());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    inline HermBandMatrix(const GenDiagMatrix<RealType(T)>& m2) :
      NEW_SIZE(m2.size(),0)
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      Zero();
      m2.AssignToD(DiagMatrixViewOf(diag()));
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

    inline HermBandMatrix(const GenDiagMatrix<ComplexType(T)>& m2) :
      NEW_SIZE(m2.size(),0)
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(IsComplex(T()));
      TMVAssert(S==RowMajor || S==ColMajor || S==DiagMajor);
      Zero();
      m2.AssignToD(DiagMatrixViewOf(diag()));
      diag().Imag().Zero();
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
    }

#undef NEW_SIZE

    virtual inline ~HermBandMatrix() 
    {
#ifdef TMVDEBUG
      SetAllTo(T(999));
#endif
    }

    //
    // Op=
    //

    inline HermBandMatrix<T,U,S,I>& operator=(
        const HermBandMatrix<T,U,S,I>& m2)
    { 
      TMVAssert(size() == m2.size());
      TMVAssert(nlo() >= m2.nlo());
      if (&m2 != this) {
        if (nlo() > m2.nlo()) LowerBand().Diags(-nlo(),-m2.nlo()).Zero();
        if (S==DiagMajor)
          if (nlo() > m2.nlo())
            std::copy(m2.start_mem(),m2.start_mem()+m2.mem_used(),
                itsm+m2.nlo()*stepi());
          else
            std::copy(m2.start_mem(),m2.start_mem()+m2.mem_used(),
                itsm1.get());
        else if (nlo()==m2.nlo())
          std::copy(m2.cptr(),m2.cptr()+m2.mem_used(),itsm1.get());
        else
          LowerBand().Diags(-m2.nlo(),1) = m2.LowerBand();
      }
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    template <IndexStyle I2> 
    inline HermBandMatrix<T,U,S,I>& operator=(
        const HermBandMatrix<T,U,S,I2>& m2)
    { 
      TMVAssert(size() == m2.size());
      TMVAssert(nlo() >= m2.nlo());
      if (nlo() > m2.nlo()) LowerBand().Diags(-nlo(),-m2.nlo()).Zero();
      if (S==DiagMajor)
        if (nlo() > m2.nlo())
          std::copy(m2.start_mem(),m2.start_mem()+m2.mem_used(),
              itsm+m2.nlo()*stepi());
        else
          std::copy(m2.start_mem(),m2.start_mem()+m2.mem_used(),
              itsm1.get());
      else if (nlo()==m2.nlo())
        std::copy(m2.cptr(),m2.cptr()+m2.mem_used(),itsm1.get());
      else
        LowerBand().Diags(-m2.nlo(),1) = m2.LowerBand();
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    inline HermBandMatrix<T,U,S,I>& operator=(
        const GenSymBandMatrix<RealType(T)>& m2)
    { 
      TMVAssert(size() == m2.size());
      TMVAssert(nlo() >= m2.nlo());
      TMVAssert(m2.isherm());
      m2.AssignTosB(View());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this; 
    }

    inline HermBandMatrix<T,U,S,I>& operator=(
        const GenSymBandMatrix<ComplexType(T)>& m2)
    { 
      TMVAssert(IsComplex(T()));
      TMVAssert(size() == m2.size());
      TMVAssert(nlo() >= m2.nlo());
      TMVAssert(m2.isherm());
      m2.AssignTosB(View());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this; 
    }

    template <class T2> 
    inline HermBandMatrix<T,U,S,I>& operator=(
        const GenSymBandMatrix<T2>& m2)
    { 
      TMVAssert(size() == m2.size());
      TMVAssert(IsReal(T2()) || m2.isherm());
      LowerBand() = m2.LowerBand();
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this; 
    }

    inline HermBandMatrix<T,U,S,I>& operator=(T x) 
    { return SetToIdentity(x); }

    inline HermBandMatrix<T,U,S,I>& operator=(
        const AssignableToSymBandMatrix<RealType(T)>& m2)
    { 
      TMVAssert(size() == m2.size());
      TMVAssert(m2.isherm());
      m2.AssignTosB(View());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    inline HermBandMatrix<T,U,S,I>& operator=(
        const AssignableToSymBandMatrix<ComplexType(T)>& m2)
    { 
      TMVAssert(IsComplex(T()));
      TMVAssert(size() == m2.size());
      TMVAssert(m2.isherm());
      m2.AssignTosB(View());
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    inline HermBandMatrix<T,U,S,I>& operator=(
        const GenDiagMatrix<RealType(T)>& m2) 
    { 
      TMVAssert(size() == m2.size());
      View() = m2;
#ifdef XTEST
      TMVAssert(this->HermOK());
#endif
      return *this;
    }

    inline HermBandMatrix<T,U,S,I>& operator=(
        const GenDiagMatrix<ComplexType(T)>& m2) 
    { TMVAssert(FALSE); return *this; }

    //
    // Access
    //

    typedef RefType(T) reference;

    inline T operator()(int i, int j) const
    { 
      if (I==CStyle) {
        TMVAssert(i>=0 && i<int(size()));
        TMVAssert(j>=0 && j<int(size()));
        return okij(i,j) ? cref(i,j) : T(0);
      } else {
        TMVAssert(i>0 && i<=int(size()));
        TMVAssert(j>0 && j<=int(size()));
        return okij(i-1,j-1) ? cref(i-1,j-1) : T(0);
      }
    }

    inline reference operator()(int i, int j) 
    { 
      if (I==CStyle) {
        TMVAssert(i>=0 && i<int(size()));
        TMVAssert(j>=0 && j<int(size()));
        TMVAssert(okij(i,j));
        return ref(i,j); 
      } else {
        TMVAssert(i>0 && i<=int(size()));
        TMVAssert(j>0 && j<=int(size()));
        TMVAssert(okij(i-1,j-1));
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
      TMVAssert(j1==j2 || okij(i,j1));
      TMVAssert(j1==j2 || okij(i,j2-1));
      TMVAssert( i-j1<=0 || j2-i<=1 );
      if ((U==Upper && i-j1<=0) || (U==Lower && j2-i<=1))
        return ConstVectorView<T,I>(itsm+i*stepi()+j1*stepj(),
            j2-j1,stepj(),NonConj); 
      else
        return ConstVectorView<T,I>(itsm+i*stepj()+j1*stepi(),
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
      TMVAssert(i1==i2 || okij(i1,j));
      TMVAssert(i1==i2 || okij(i2-1,j));
      TMVAssert( j-i1<=0 || i2-j<=1 );
      if ((U==Upper && i2-j<=1) || (U==Lower && j-i1<=0))
        return ConstVectorView<T,I>(itsm+i1*stepi()+j*stepj(),
            i2-i1,stepi(),NonConj); 
      else 
        return ConstVectorView<T,I>(itsm+i1*stepj()+j*stepi(),
            i2-i1,stepj(),ConjOf(T,NonConj)); 
    }

    inline ConstVectorView<T,I> diag() const
    { return ConstVectorView<T,I>(itsm,size(),diagstep(),NonConj); }

    inline ConstVectorView<T,I> diag(int i) const
    {
      TMVAssert(i>=-nlo() && i<=nlo());
      TMVAssert(i>=-int(size()) && i<=int(size())); 
      ConjType newct = 
      ((i>0) == (U==Upper)) ? NonConj : ConjOf(T,NonConj); 
      i = std::abs(i);
      return ConstVectorView<T,I>(itsm+i*(U==Upper?stepj():stepi()),
          size()-i,diagstep(),newct); 
    }

    inline ConstVectorView<T,I> diag(int i, int j1, int j2) const
    {
      TMVAssert(i>=-nlo() && i<=nlo());
      TMVAssert(i>=-int(size()) && i<=int(size())); 
      if (I==FortranStyle) { 
        TMVAssert(j1>0 && j1-j2<=0 && j2<=int(size())-std::abs(i)); --j1;
      } else {
        TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size())-std::abs(i));
      }
      ConjType newct = 
      ((i>0) == (U==Upper)) ? NonConj : ConjOf(T,NonConj); 
      i = std::abs(i);
      return ConstVectorView<T,I>(
          itsm+i*(U==Upper?stepj():stepi())+j1*diagstep(),
          j2-j1,diagstep(),newct);
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
      TMVAssert(j1==j2 || okij(i,j1));
      TMVAssert(j1==j2 || okij(i,j2-1));
      TMVAssert( i-j1<=0 || j2-i<=1 );
      if ((U==Upper && i-j1<=0) || (U==Lower && j2-i<=1))
        return VectorView<T,I>(itsm+i*stepi()+j1*stepj(),
            j2-j1,stepj(),NonConj FIRSTLAST); 
      else
        return VectorView<T,I>(itsm+i*stepj()+j1*stepi(),
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
      TMVAssert(i1==i2 || okij(i1,j));
      TMVAssert(i1==i2 || okij(i2-1,j));
      TMVAssert( j-i1<=0 || i2-j<=1 );
      if ((U==Upper && i2-j<=1) || (U==Lower && j-i1<=0))
        return VectorView<T,I>(itsm+i1*stepi()+j*stepj(),
            i2-i1,stepi(),NonConj FIRSTLAST); 
      else 
        return VectorView<T,I>(itsm+i1*stepj()+j*stepi(),
            i2-i1,stepj(),ConjOf(T,NonConj) FIRSTLAST); 
    }

    inline VectorView<T,I> diag()
    { return VectorView<T,I>(itsm,size(),diagstep(),NonConj FIRSTLAST); }

    inline VectorView<T,I> diag(int i) 
    {
      TMVAssert(i>=-nlo() && i<=nlo());
      TMVAssert(i>=-int(size()) && i<=int(size())); 
      ConjType newct = 
      ((i>0) == (U==Upper)) ? NonConj : ConjOf(T,NonConj); 
      i = std::abs(i);
      return VectorView<T,I>(itsm+i*(U==Upper?stepj():stepi()),
          size()-i,diagstep(),newct FIRSTLAST); 
    }

    inline VectorView<T,I> diag(int i, int j1, int j2) 
    {
      TMVAssert(i>=-nlo() && i<=nlo());
      TMVAssert(i>=-int(size()) && i<=int(size())); 
      if (I==FortranStyle) { 
        TMVAssert(j1>0 && j1-j2<=0 && j2<=int(size())-std::abs(i)); --j1;
      } else {
        TMVAssert(j1>=0 && j1-j2<=0 && j2<=int(size())-std::abs(i));
      }
      ConjType newct = 
      ((i>0) == (U==Upper)) ? NonConj : ConjOf(T,NonConj); 
      i = std::abs(i);
      return VectorView<T,I>(
          itsm+i*(U==Upper?stepj():stepi())+j1*diagstep(),
          j2-j1,diagstep(),newct FIRSTLAST);
    }


    //
    // Modifying Functions
    //

    inline HermBandMatrix<T,U,S,I>& Zero() 
    {
      std::fill_n(itsm1.get(),linsize,T(0));
      return *this; 
    }

    inline HermBandMatrix<T,U,S,I>& SetAllTo(T x) 
    {
      TMVAssert(IMAG(x) == RealType(T)(0));
      UpperBand().SetAllTo(x); 
      return *this; 
    }

    inline HermBandMatrix<T,U,S,I>& Clip(RealType(T) thresh) 
    { UpperBand().Clip(thresh); return *this; }

    inline HermBandMatrix<T,U,S,I>& ConjugateSelf() 
    { if (IsComplex(T())) UpperBandOff().ConjugateSelf(); return *this; }

    inline HermBandMatrix<T,U,S,I>& TransposeSelf() 
    { if (IsComplex(T())) UpperBandOff().ConjugateSelf(); return *this; }

    inline HermBandMatrix<T,U,S,I>& SetToIdentity(T x=T(1)) 
    {
      TMVAssert(IMAG(x) == RealType(T)(0));
      Zero(); 
      diag().SetAllTo(x); 
      return *this; 
    }

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
        return ConstMatrixView<T,I>(itsm+i1*stepi()+j1*stepj(),
            i2-i1,j2-j1,stepi(),stepj(),S,NonConj);
      else
        return ConstMatrixView<T,I>(itsm+i1*stepj()+j1*stepi(),
            i2-i1,j2-j1,stepj(),stepi(),TransOf(S),ConjOf(T,NonConj));
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
        return ConstMatrixView<T,I>(itsm+i1*stepi()+j1*stepj(),
            (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
            newstor,NonConj);
      else
        return ConstMatrixView<T,I>(itsm+i1*stepj()+j1*stepi(),
            (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
            TransOf(newstor),ConjOf(T,NonConj));
    }

    inline ConstBandMatrixView<T,I> SubBandMatrix(
        int i1, int i2, int j1, int j2, int newnlo, int newnhi) const
    {
      TMVAssert(View().OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
      if (I==FortranStyle) { --i1; --j1; }
      TMVAssert(i1+newnlo-j1<=0 || j1+newnhi-i1<=0);
      if ((i1+newnlo-j1<=0 && U==Upper) || (j1+newnhi-i1<=0 && U==Lower))
        return ConstBandMatrixView<T,I>(itsm+i1*stepi()+j1*stepj(),
            i2-i1,j2-j1,newnlo,newnhi,stepi(),stepj(),diagstep(),
            S,NonConj);
      else
        return ConstBandMatrixView<T,I>(itsm+i1*stepj()+j1*stepi(),
            i2-i1,j2-j1,newnlo,newnhi,stepj(),stepi(),diagstep(),
            TransOf(S),ConjOf(T,NonConj));
    }

    inline ConstBandMatrixView<T,I> SubBandMatrix(
        int i1, int i2, int j1, int j2, int newnlo, int newnhi,
        int istep, int jstep) const
    {
      TMVAssert(View().OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,
            istep,jstep));
      if (I==FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
      const StorageType newstor=
      iscm() ? (istep == 1 ? ColMajor : NoMajor) :
      isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
      TMVAssert(i1+newnlo*istep<=j1 || j1+newnhi*jstep<=i1);
      if ((i1+newnlo*istep<=j1 && U==Upper) || 
          (j1+newnhi*jstep<=i1 && U==Lower)) {
        const int newstepi = stepi()*istep;
        const int newstepj = stepj()*jstep;
        return ConstBandMatrixView<T,I>(itsm+i1*stepi()+j1*stepj(),
            (i2-i1)/istep,(j2-j1)/jstep,newnlo,newnhi,
            newstepi,newstepj,newstepi+newstepj,newstor,NonConj);
      } else {
        const int newstepi = stepj()*istep;
        const int newstepj = stepi()*jstep;
        return ConstBandMatrixView<T,I>(itsm+i1*stepj()+j1*stepi(),
            (i2-i1)/istep,(j2-j1)/jstep,newnlo,newnhi,
            newstepi,newstepj,newstepi+newstepj,
            TransOf(newstor),ConjOf(T,NonConj));
      }
    }

    inline ConstVectorView<T,I> SubVector(int i, int j,
        int istep, int jstep, int n) const
    {
      TMVAssert(View().OKSubVector(i,j,istep,jstep,n));
      if (I==FortranStyle) { --i; --j; }
      if ((U==Upper && i-j<=0) || (U==Lower && j-i<=0))
        return ConstVectorView<T,I>(itsm+i*stepi()+j*stepj(),n,
            istep*stepi()+jstep*stepj(),NonConj);
      else
        return ConstVectorView<T,I>(itsm+i*stepj()+j*stepi(),n,
            istep*stepj()+jstep*stepi(),ConjOf(T,NonConj));
    }

    inline ConstSymMatrixView<T,I> SubSymMatrix(int i1, int i2) const
    {
      TMVAssert(View().OKSubSymMatrix(i1,i2,1));
      if (I==FortranStyle) { --i1; }
      return ConstSymMatrixView<T,I>(itsm+i1*diagstep(),i2-i1,
          stepi(),stepj(),Herm,U,S,NonConj);
    }

    inline ConstSymMatrixView<T,I> SubSymMatrix(int i1, int i2,
        int istep) const
    {
      TMVAssert(View().OKSubSymMatrix(i1,i2,istep));
      if (I==FortranStyle) { --i1; i2+=istep-1; }
      return ConstSymMatrixView<T,I>(itsm+i1*diagstep(),
          (i2-i1)/istep,istep*stepi(),istep*stepj(),Herm,U,
          istep==1 ? S : NoMajor, NonConj);
    }

    inline ConstSymBandMatrixView<T,I> SubSymBandMatrix(
        int i1, int i2, int newnlo) const
    {
      TMVAssert(View().OKSubSymBandMatrix(i1,i2,newnlo,1));
      if (I==FortranStyle) { --i1; }
      return ConstSymBandMatrixView<T,I>(itsm+i1*diagstep(),i2-i1,
          newnlo,stepi(),stepj(),diagstep(),Herm,U,S,NonConj);
    }

    inline ConstSymBandMatrixView<T,I> SubSymBandMatrix(
        int i1, int i2) const
    { 
      return SubSymBandMatrix(i1,i2,MIN(nlo(),i2-i1-(I==CStyle?0:1))); 
    }

    inline ConstSymBandMatrixView<T,I> SubSymBandMatrix(
        int i1, int i2, int newnlo, int istep) const
    {
      TMVAssert(View().OKSubSymBandMatrix(i1,i2,newnlo,istep));
      if (I==FortranStyle) { --i1; i2+=istep-1; }
      return ConstSymBandMatrixView<T,I>(
          itsm+i1*diagstep(), (i2-i1)/istep, newnlo,
          istep*stepi(),istep*stepj(),istep*diagstep(),
          Herm,U,istep==1 ? S : NoMajor,NonConj);
    }

    inline ConstSymBandMatrixView<T,I> SymDiags(int newnlo) const
    {
      TMVAssert(newnlo>=0 && newnlo <= nlo());
      return ConstSymBandMatrixView<T,I>(itsm,size(),newnlo,
          stepi(),stepj(),diagstep(),Herm,U,S,NonConj);
    }

    inline ConstBandMatrixView<T> Diags(int k1, int k2) const
    {
      TMVAssert(k1>=-nlo() && k1<k2 && k2<=nlo()+1);
      TMVAssert(k2 <= 1 || k1 >= 0);

      if (k1 >= 0) {
        const int newsize = int(size())-k1;
        const int newnhi = k2-k1-1;
        if (uplo() == Upper)
          return ConstBandMatrixView<T>(itsm+k1*stepj(),
              newsize, newsize, 0, newnhi, stepi(), stepj(), diagstep(),
              stor(), ct());
        else
          return ConstBandMatrixView<T>(itsm+k1*stepi(),
              newsize, newsize, 0, newnhi, stepj(), stepi(), diagstep(),
              TransOf(stor()), issym()?ct():ConjOf(T,ct()));
      } else {
        const int newsize = int(size())+k2-1;
        const int newnlo = k2-k1-1;
        if (uplo() == Lower)
          return ConstBandMatrixView<T>(itsm-k2*stepi(),
              newsize, newsize, newnlo, 0, stepi(), stepj(), diagstep(),
              stor(), ct());
        else
          return ConstBandMatrixView<T>(itsm-k2*stepj(),
              newsize, newsize, newnlo, 0, stepj(), stepi(), diagstep(),
              TransOf(stor()), issym()?ct():ConjOf(T,ct()));
      }
    }

    inline ConstBandMatrixView<T> UpperBand() const
    {
      if (uplo() == Upper)
        return ConstBandMatrixView<T>(itsm,size(),size(),0,nlo(),
            stepi(),stepj(),diagstep(),S,NonConj);
      else
        return ConstBandMatrixView<T>(itsm,size(),size(),0,nlo(),
            stepj(),stepi(),diagstep(),TransOf(S),ConjOf(T,NonConj));
    }

    inline ConstBandMatrixView<T> UpperBandOff() const
    {
      if (uplo() == Upper)
        return ConstBandMatrixView<T>(itsm+stepj(),size()-1,size()-1,
            0,nlo()-1,stepi(),stepj(),diagstep(),S,NonConj);
      else
        return ConstBandMatrixView<T>(itsm+stepi(),size()-1,size()-1,
            0,nlo()-1,stepj(),stepi(),diagstep(),
            TransOf(S),ConjOf(T,NonConj));
    }

    inline ConstBandMatrixView<T> LowerBand() const
    {
      if (uplo() == Lower)
        return ConstBandMatrixView<T>(itsm,size(),size(),nlo(),0,
            stepi(),stepj(),diagstep(),S,NonConj);
      else
        return ConstBandMatrixView<T>(itsm,size(),size(),nlo(),0,
            stepj(),stepi(),diagstep(),TransOf(S),ConjOf(T,NonConj));
    }

    inline ConstBandMatrixView<T> LowerBandOff() const
    {
      if (uplo() == Lower)
        return ConstBandMatrixView<T>(itsm+stepi(),size()-1,size()-1,
            nlo()-1,0,stepi(),stepj(),diagstep(),S,NonConj);
      else
        return ConstBandMatrixView<T>(itsm+stepj(),size()-1,size()-1,
            nlo()-1,0,stepj(),stepi(),diagstep(),
            TransOf(S),ConjOf(T,NonConj));
    }

    inline MatrixView<T,I> SubMatrix(int i1, int i2, int j1, int j2)
    {
      TMVAssert(View().OKSubMatrix(i1,i2,j1,j2,1,1));
      if (I==FortranStyle) { --i1; --j1; }
      TMVAssert(i2-1<=j1 || j2-1<=i1);
      if ((U==Upper && i2-1<=j1) || (U==Lower && j2-1<=i1))
        return MatrixView<T,I>(itsm+i1*stepi()+j1*stepj(),
            i2-i1,j2-j1,stepi(),stepj(),S,NonConj FIRSTLAST);
      else
        return MatrixView<T,I>(itsm+i1*stepj()+j1*stepi(),
            i2-i1,j2-j1,stepj(),stepi(),TransOf(S),ConjOf(T,NonConj) 
            FIRSTLAST);
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
        return MatrixView<T,I>(itsm+i1*stepi()+j1*stepj(),
            (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
            newstor,NonConj FIRSTLAST);
      else
        return MatrixView<T,I>(itsm+i1*stepj()+j1*stepi(),
            (i2-i1)/istep,(j2-j1)/jstep,istep*stepj(),jstep*stepi(),
            TransOf(newstor),ConjOf(T,NonConj) FIRSTLAST);
    }

    inline BandMatrixView<T,I> SubBandMatrix(
        int i1, int i2, int j1, int j2, int newnlo, int newnhi)
    {
      TMVAssert(View().OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,1,1));
      if (I==FortranStyle) { --i1; --j1; }
      TMVAssert(i1+newnlo-j1<=0 || j1+newnhi-i1<=0);
      if ((i1+newnlo-j1<=0 && U==Upper) || (j1+newnhi-i1<=0 && U==Lower))
        return BandMatrixView<T,I>(itsm+i1*stepi()+j1*stepj(),
            i2-i1,j2-j1,newnlo,newnhi,stepi(),stepj(),diagstep(),
            S,NonConj FIRSTLAST);
      else
        return BandMatrixView<T,I>(itsm+i1*stepj()+j1*stepi(),
            i2-i1,j2-j1,newnlo,newnhi,stepj(),stepi(),diagstep(),
            TransOf(S),ConjOf(T,NonConj) FIRSTLAST);
    }

    inline BandMatrixView<T,I> SubBandMatrix(
        int i1, int i2, int j1, int j2, int newnlo, int newnhi,
        int istep, int jstep)
    {
      TMVAssert(View().OKSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,
            istep,jstep));
      if (I==FortranStyle) { --i1; --j1; i2+=istep-1; j2+=jstep-1; }
      const StorageType newstor=
      iscm() ? (istep == 1 ? ColMajor : NoMajor) :
      isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
      TMVAssert(i1+newnlo*istep<=j1 || j1+newnhi*jstep<=i1);
      if ((i1+newnlo*istep<=j1 && U==Upper) || 
          (j1+newnhi*jstep<=i1 && U==Lower)) {
        const int newstepi = stepi()*istep;
        const int newstepj = stepj()*jstep;
        return BandMatrixView<T,I>(itsm+i1*stepi()+j1*stepj(),
            (i2-i1)/istep,(j2-j1)/jstep,newnlo,newnhi,
            newstepi,newstepj,newstepi+newstepj,newstor,NonConj FIRSTLAST);
      } else {
        const int newstepi = stepj()*istep;
        const int newstepj = stepi()*jstep;
        return BandMatrixView<T,I>(itsm+i1*stepj()+j1*stepi(),
            (i2-i1)/istep,(j2-j1)/jstep,newnlo,newnhi,
            newstepi,newstepj,newstepi+newstepj,
            TransOf(newstor),ConjOf(T,NonConj) FIRSTLAST);
      }
    }

    inline VectorView<T,I> SubVector(int i, int j,
        int istep, int jstep, int n)
    {
      TMVAssert(View().OKSubVector(i,j,istep,jstep,n));
      if (I==FortranStyle) { --i; --j; }
      if ((U==Upper && i-j<=0) || (U==Lower && j-i<=0))
        return VectorView<T,I>(itsm+i*stepi()+j*stepj(),n,
            istep*stepi()+jstep*stepj(),NonConj FIRSTLAST);
      else
        return VectorView<T,I>(itsm+i*stepj()+j*stepi(),n,
            istep*stepj()+jstep*stepi(),ConjOf(T,NonConj) FIRSTLAST);
    }

    inline SymMatrixView<T,I> SubSymMatrix(int i1, int i2)
    {
      TMVAssert(View().OKSubSymMatrix(i1,i2,1));
      if (I==FortranStyle) { --i1; }
      return SymMatrixView<T,I>(itsm+i1*diagstep(),i2-i1,
          stepi(),stepj(),Herm,U,S,NonConj FIRSTLAST);
    }

    inline SymMatrixView<T,I> SubSymMatrix(int i1, int i2, int istep) 
    {
      TMVAssert(View().OKSubSymMatrix(i1,i2,istep));
      if (I==FortranStyle) { --i1; i2+=istep-1; }
      return SymMatrixView<T,I>(itsm+i1*diagstep(),(i2-i1)/istep,
          istep*stepi(),istep*stepj(),Herm,U,
          istep==1 ? S : NoMajor,NonConj FIRSTLAST);
    }

    inline SymBandMatrixView<T,I> SubSymBandMatrix(int i1, int i2, int newnlo)
    {
      TMVAssert(View().OKSubSymBandMatrix(i1,i2,newnlo,1));
      if (I==FortranStyle) { --i1; }
      return SymBandMatrixView<T,I>(itsm+i1*diagstep(),i2-i1,newnlo,
          stepi(),stepj(),diagstep(),Herm,U,S,NonConj FIRSTLAST);
    }

    inline SymBandMatrixView<T,I> SubSymBandMatrix(int i1, int i2)
    { 
      return SubSymBandMatrix(i1,i2,MIN(nlo(),i2-i1-(I==CStyle?0:1))); 
    }

    inline SymBandMatrixView<T,I> SubSymBandMatrix(
        int i1, int i2, int newnlo, int istep)
    {
      TMVAssert(View().OKSubSymBandMatrix(i1,i2,newnlo,istep));
      if (I==FortranStyle) { --i1; i2+=istep-1; }
      return SymBandMatrixView<T,I>(itsm+i1*diagstep(),(i2-i1)/istep,
          newnlo,istep*stepi(),istep*stepj(),istep*diagstep(),
          Herm,U,istep==1 ? S : NoMajor,NonConj FIRSTLAST);
    }

    inline SymBandMatrixView<T,I> SymDiags(int newnlo)
    {
      TMVAssert(newnlo>=0 && newnlo <= nlo());
      return SymBandMatrixView<T,I>(itsm,size(),newnlo,
          stepi(),stepj(),diagstep(),Herm,U,S,NonConj FIRSTLAST);
    }

    inline BandMatrixView<T> Diags(int k1, int k2)
    {
      TMVAssert(k1>=-nlo() && k1<k2 && k2<=nlo()+1);
      TMVAssert(k2 <= 1 || k1 >= 0);

      if (k1 >= 0) {
        const int newsize = int(size())-k1;
        const int newnhi = k2-k1-1;
        if (uplo() == Upper)
          return BandMatrixView<T>(itsm+k1*stepj(),
              newsize, newsize, 0, newnhi, stepi(), stepj(), diagstep(),
              stor(), ct() FIRSTLAST );
        else
          return BandMatrixView<T>(itsm+k1*stepi(),
              newsize, newsize, 0, newnhi, stepj(), stepi(), diagstep(),
              TransOf(stor()), issym()?ct():ConjOf(T,ct()) FIRSTLAST );
      } else {
        const int newsize = int(size())+k2-1;
        const int newnlo = k2-k1-1;
        if (uplo() == Lower)
          return BandMatrixView<T>(itsm-k2*stepi(),
              newsize, newsize, newnlo, 0, stepi(), stepj(), diagstep(),
              stor(), ct() FIRSTLAST );
        else
          return BandMatrixView<T>(itsm-k2*stepj(),
              newsize, newsize, newnlo, 0, stepj(), stepi(), diagstep(),
              TransOf(stor()), issym()?ct():ConjOf(T,ct()) FIRSTLAST );
      }
    }

    inline BandMatrixView<T> UpperBand()
    {
      if (uplo() == Upper)
        return BandMatrixView<T>(itsm,size(),size(),0,nlo(),
            stepi(),stepj(),diagstep(),S,NonConj FIRSTLAST);
      else
        return BandMatrixView<T>(itsm,size(),size(),0,nlo(),
            stepj(),stepi(),diagstep(),TransOf(S),ConjOf(T,NonConj) 
            FIRSTLAST);
    }

    inline BandMatrixView<T> UpperBandOff()
    {
      if (uplo() == Upper)
        return BandMatrixView<T>(itsm+stepj(),size()-1,size()-1,
            0,nlo()-1,stepi(),stepj(),diagstep(),S,NonConj FIRSTLAST);
      else
        return BandMatrixView<T>(itsm+stepi(),size()-1,size()-1,
            0,nlo()-1,stepj(),stepi(),diagstep(),
            TransOf(S),ConjOf(T,NonConj) FIRSTLAST);
    }

    inline BandMatrixView<T> LowerBand()
    {
      if (uplo() == Lower)
        return BandMatrixView<T>(itsm,size(),size(),nlo(),0,
            stepi(),stepj(),diagstep(),S,NonConj FIRSTLAST);
      else
        return BandMatrixView<T>(itsm,size(),size(),nlo(),0,
            stepj(),stepi(),diagstep(),TransOf(S),ConjOf(T,NonConj) 
            FIRSTLAST);
    }

    inline BandMatrixView<T> LowerBandOff()
    {
      if (uplo() == Lower)
        return BandMatrixView<T>(itsm+stepi(),size()-1,size()-1,
            nlo()-1,0,stepi(),stepj(),diagstep(),S,NonConj FIRSTLAST);
      else
        return BandMatrixView<T>(itsm+stepj(),size()-1,size()-1,
            nlo()-1,0,stepj(),stepi(),diagstep(),
            TransOf(S),ConjOf(T,NonConj) FIRSTLAST);
    }

    inline ConstSymBandMatrixView<T,I> View() const
    { 
      return ConstSymBandMatrixView<T,I>(itsm,size(),nlo(),
          stepi(),stepj(),diagstep(),Herm,U,S,NonConj);
    }

    inline ConstSymBandMatrixView<T,I> Transpose() const
    {
      return ConstSymBandMatrixView<T,I>(itsm,size(),nlo(),
          stepj(),stepi(),diagstep(),Herm,UTransOf(U),TransOf(S),NonConj);
    }

    inline ConstSymBandMatrixView<T,I> Conjugate() const
    { 
      return ConstSymBandMatrixView<T,I>(itsm,size(),nlo(),
          stepi(),stepj(),diagstep(),Herm,U,S,ConjOf(T,NonConj));
    }

    inline ConstSymBandMatrixView<T,I> Adjoint() const
    {
      return ConstSymBandMatrixView<T,I>(itsm,size(),nlo(),
          stepj(),stepi(),diagstep(),
          Herm,UTransOf(U),TransOf(S),ConjOf(T,NonConj));
    }

    inline ConstSymBandMatrixView<RealType(T),I> Real() const
    {
      return ConstSymBandMatrixView<RealType(T)>(
          reinterpret_cast<const RealType(T)*>(itsm),size(),nlo(),
          IsReal(T()) ? stepi() : 2*stepi(),
          IsReal(T()) ? stepj() : 2*stepj(),
          IsReal(T()) ? diagstep() : 2*diagstep(),
          Herm,U,IsReal(T())?S:NoMajor,NonConj);
    }

    inline ConstSymBandMatrixView<RealType(T),I> Imag() const
    {
      TMVAssert(IsComplex(T()));
      return ConstSymBandMatrixView<RealType(T)>(
          reinterpret_cast<const RealType(T)*>(itsm)+1,size(),nlo(),
          2*stepi(),2*stepj(),2*diagstep(),Herm,U,NoMajor,NonConj);
    } 

    inline SymBandMatrixView<T,I> View() 
    { 
      return SymBandMatrixView<T,I>(itsm,size(),nlo(),
          stepi(),stepj(),diagstep(), Herm,U,S,NonConj FIRSTLAST);
    }

    inline SymBandMatrixView<T,I> Transpose() 
    {
      return SymBandMatrixView<T,I>(itsm,size(),nlo(),
          stepj(),stepi(),diagstep(),
          Herm,UTransOf(U),TransOf(S),NonConj FIRSTLAST);
    }

    inline SymBandMatrixView<T,I> Conjugate() 
    { 
      return SymBandMatrixView<T,I>(itsm,size(),nlo(),
          stepi(),stepj(),diagstep(),
          Herm,U,S,ConjOf(T,NonConj) FIRSTLAST);
    }

    inline SymBandMatrixView<T,I> Adjoint() 
    {
      return SymBandMatrixView<T,I>(itsm,size(),nlo(),
          stepj(),stepi(),diagstep(),
          Herm,UTransOf(U),TransOf(S),ConjOf(T,NonConj) FIRSTLAST);
    }

    inline SymBandMatrixView<RealType(T),I> Real()
    {
      return SymBandMatrixView<RealType(T)>(
          reinterpret_cast<RealType(T)*>(itsm),size(),nlo(),
          IsReal(T()) ? stepi() : 2*stepi(),
          IsReal(T()) ? stepj() : 2*stepj(),
          IsReal(T()) ? diagstep() : 2*diagstep(),
          Herm,U,IsReal(T())?S:NoMajor,NonConj
#ifdef TMVFLDEBUG
          ,reinterpret_cast<const RealType(T)*>(first)
          ,reinterpret_cast<const RealType(T)*>(last)
#endif
          );
    }

    inline SymBandMatrixView<RealType(T),I> Imag()
    // The imaginary part of a Hermitian matrix is anti-symmetric
    // so this is illegal.
    {
      TMVAssert(FALSE);
      return SymBandMatrixView<RealType(T),I>(0,0,0,0,0,0,Herm,U,S,NonConj
          FIRSTLAST1(0,0) );
    } 

    inline size_t size() const { return itss; }
    inline int nlo() const { return itslo; }
    inline size_t mem_used() const { return linsize; }
    inline const T* start_mem() const { return itsm1.get(); }
    inline const T* cptr() const { return itsm; }
    inline T* ptr() { return itsm; }
    inline int stepi() const { return itssi; }
    inline int stepj() const { return itssj; }
    inline int diagstep() const { return itssd; }
    inline SymType sym() const { return Herm; }
    inline UpLoType uplo() const { return U; }
    inline StorageType stor() const { return S; }
    inline ConjType ct() const { return NonConj; }
    inline bool isrm() const { return S==RowMajor; }
    inline bool iscm() const { return S==ColMajor; }
    inline bool isdm() const { return S==DiagMajor; }
    inline bool isconj() const { return false; }
    inline bool isherm() const { return true; }
    inline bool issym() const { return IsReal(T()); }

    inline reference ref(int i, int j)
    {
      if ((U==Upper && i <= j) || (U==Lower && i>=j)) 
        return REF(itsm + i*stepi() + j*stepj(),NonConj);
      else 
        return REF(itsm + j*stepi() + i*stepj(),Conj);
    }

    inline T cref(int i, int j) const 
    {
      if ((U==Upper && i <= j) || (U==Lower && i>=j)) 
        return itsm[i*itssi + j*itssj];
      else 
        return CONJ(itsm[j*itssi + i*itssj]);
    }

  protected :

    const size_t linsize;
    auto_array<T> itsm1;
    const size_t itss;
    const int itslo;
    const int itssi;
    const int itssj;
    const int itssd;
    T*const itsm;

#ifdef TMVFLDEBUG
  public:
    const T*const first;
    const T*const last;
  protected:
#endif

    inline bool okij(int i, int j) const
    { return (j+nlo() >= i && i+nlo() >= j); }

  }; // HermBandMatrix

  //---------------------------------------------------------------------------

  //
  // Special Creators: 
  //   SymTriDiagMatrix(v1,v2)
  //       Returns a SymBandMatrix with v1 as the diagonal, and v2 as
  //       the offdiagonal.
  //   HermTriDiagMatrix(v1,v2,uplo)
  //       Returns a HermBandMatrix with v1 as the diagonal, and v2 as
  //       the upper or lower offdiagonal, according to the value of uplo.
  //   SymBandMatrixViewOf(m,uplo,nlo)
  //   HermBandMatrixViewOf(m,uplo,nlo)
  //   SymBandMatrixViewOf(mptr,size,nlo,uplo,stor)
  //   HermBandMatrixViewOf(mptr,size,nlo,uplo,stor)
  //

  template <class T> 
  SymBandMatrix<T,Upper,DiagMajor> SymTriDiagMatrix(
      const GenVector<T>& v1, const GenVector<T>& v2);
  template <class T> 
  HermBandMatrix<T,Upper,DiagMajor> HermTriDiagMatrix(
      const GenVector<T>& v1, const GenVector<T>& v2, 
      UpLoType uplo=tmv::Upper);
  template <class T> 
  HermBandMatrix<std::complex<T>,Upper,DiagMajor> HermTriDiagMatrix(
      const GenVector<T>& v1, const GenVector<std::complex<T> >& v2,
      UpLoType uplo);

  // From Matrix
  template <class T> 
  inline ConstSymBandMatrixView<T> SymBandMatrixViewOf(
      const GenMatrix<T>& m, UpLoType uplo, int nlo)
  {
    TMVAssert(m.colsize()==m.rowsize());
    return ConstSymBandMatrixView<T>(m.cptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),Sym,uplo,m.stor(),m.ct()); 
  }

  template <class T, IndexStyle I> 
  inline ConstSymBandMatrixView<T,I> SymBandMatrixViewOf(
      const ConstMatrixView<T,I>& m, UpLoType uplo, int nlo)
  { 
    TMVAssert(m.colsize()==m.rowsize());
    return ConstSymBandMatrixView<T,I>(m.cptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),Sym,uplo,m.stor(),m.ct()); 
  }

  template <class T, StorageType S, IndexStyle I> 
  inline ConstSymBandMatrixView<T,I> SymBandMatrixViewOf(
      const Matrix<T,S,I>& m, UpLoType uplo, int nlo)
  {
    TMVAssert(m.colsize()==m.rowsize());
    return ConstSymBandMatrixView<T,I>(m.cptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),Sym,uplo,m.stor(),m.ct()); 
  }

  template <class T, IndexStyle I> 
  inline SymBandMatrixView<T,I> SymBandMatrixViewOf(
      const MatrixView<T,I>& m, UpLoType uplo, int nlo)
  { 
    TMVAssert(m.colsize()==m.rowsize());
    return SymBandMatrixView<T,I>(m.ptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),Sym,uplo,m.stor(),m.ct()
        FIRSTLAST1(m.first,m.last) ); 
  }

  template <class T, StorageType S, IndexStyle I> 
  inline SymBandMatrixView<T,I> SymBandMatrixViewOf(
      Matrix<T,S,I>& m, UpLoType uplo, int nlo)
  {
    TMVAssert(m.colsize()==m.rowsize());
    return SymBandMatrixView<T,I>(m.ptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),Sym,uplo,m.stor(),m.ct()
        FIRSTLAST1(m.first,m.last) ); 
  }

  template <class T> 
  inline ConstSymBandMatrixView<T> HermBandMatrixViewOf(
      const GenMatrix<T>& m, UpLoType uplo, int nlo)
  {
    TMVAssert(m.colsize()==m.rowsize());
#ifdef XTEST_DEBUG
    TMVAssert(IsReal(T()) || m.diag().Imag().NormInf() == RealType(T)(0));
#endif
    return ConstSymBandMatrixView<T>(m.cptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),Herm,uplo,m.stor(),m.ct()); 
  }

  template <class T, IndexStyle I> 
  inline ConstSymBandMatrixView<T> HermBandMatrixViewOf(
      const ConstMatrixView<T,I>& m, UpLoType uplo, int nlo)
  { 
    TMVAssert(m.colsize()==m.rowsize());
#ifdef XTEST_DEBUG
    TMVAssert(IsReal(T()) || m.diag().Imag().NormInf() == RealType(T)(0));
#endif
    return ConstSymBandMatrixView<T,I>(m.cptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),Herm,uplo,m.stor(),m.ct()); 
  }

  template <class T, StorageType S, IndexStyle I> 
  inline ConstSymBandMatrixView<T,I> HermBandMatrixViewOf(
      const Matrix<T,S,I>& m, UpLoType uplo, int nlo)
  {
    TMVAssert(m.colsize()==m.rowsize());
#ifdef XTEST_DEBUG
    TMVAssert(IsReal(T()) || m.diag().Imag().NormInf() == RealType(T)(0));
#endif
    return ConstSymBandMatrixView<T,I>(m.cptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),Herm,uplo,m.stor(),m.ct()); 
  }

  template <class T, IndexStyle I> 
  inline SymBandMatrixView<T> HermBandMatrixViewOf(
      const MatrixView<T,I>& m, UpLoType uplo, int nlo)
  { 
    TMVAssert(m.colsize()==m.rowsize());
#ifdef XTEST_DEBUG
    TMVAssert(IsReal(T()) || m.diag().Imag().NormInf() == RealType(T)(0));
#endif
    return SymBandMatrixView<T,I>(m.ptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),Herm,uplo,m.stor(),m.ct()
        FIRSTLAST1(m.first,m.last) ); 
  }

  template <class T, StorageType S, IndexStyle I> 
  inline SymBandMatrixView<T,I> HermBandMatrixViewOf(
      Matrix<T,S,I>& m, UpLoType uplo, int nlo)
  {
    TMVAssert(m.colsize()==m.rowsize());
#ifdef XTEST_DEBUG
    TMVAssert(IsReal(T()) || m.diag().Imag().NormInf() == RealType(T)(0));
#endif
    return SymBandMatrixView<T,I>(m.ptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),Herm,uplo,m.stor(),m.ct()
        FIRSTLAST1(m.first,m.last) ); 
  }

  // From BandMatrix
  template <class T> 
  inline ConstSymBandMatrixView<T> SymBandMatrixViewOf(
      const GenBandMatrix<T>& m, UpLoType uplo, int nlo=-1)
  {
    if (nlo<0) nlo = (uplo==Upper?m.nhi():m.nlo());
    TMVAssert(m.colsize()==m.rowsize());
    return ConstSymBandMatrixView<T>(m.cptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),Sym,uplo,m.stor(),m.ct()); 
  }

  template <class T, IndexStyle I> 
  inline ConstSymBandMatrixView<T,I> SymBandMatrixViewOf(
      const ConstBandMatrixView<T,I>& m, UpLoType uplo, int nlo=-1)
  { 
    if (nlo<0) nlo = (uplo==Upper?m.nhi():m.nlo());
    TMVAssert(m.colsize()==m.rowsize());
    return ConstSymBandMatrixView<T,I>(m.cptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),Sym,uplo,m.stor(),m.ct()); 
  }

  template <class T, StorageType S, IndexStyle I> 
  inline ConstSymBandMatrixView<T,I> SymBandMatrixViewOf(
      const BandMatrix<T,S,I>& m, UpLoType uplo, int nlo=-1)
  {
    if (nlo<0) nlo = (uplo==Upper?m.nhi():m.nlo());
    TMVAssert(m.colsize()==m.rowsize());
    return ConstSymBandMatrixView<T,I>(m.cptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),Sym,uplo,m.stor(),m.ct()); 
  }

  template <class T, IndexStyle I> 
  inline SymBandMatrixView<T,I> SymBandMatrixViewOf(
      const BandMatrixView<T,I>& m, UpLoType uplo, int nlo=-1)
  { 
    if (nlo<0) nlo = (uplo==Upper?m.nhi():m.nlo());
    TMVAssert(m.colsize()==m.rowsize());
    return SymBandMatrixView<T,I>(m.ptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),Sym,uplo,m.stor(),m.ct()
        FIRSTLAST1(m.first,m.last) ); 
  }

  template <class T, StorageType S, IndexStyle I> 
  inline SymBandMatrixView<T,I> SymBandMatrixViewOf(
      BandMatrix<T,S,I>& m, UpLoType uplo, int nlo=-1)
  {
    if (nlo<0) nlo = (uplo==Upper?m.nhi():m.nlo());
    TMVAssert(m.colsize()==m.rowsize());
    return SymBandMatrixView<T,I>(m.ptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),Sym,uplo,m.stor(),m.ct()
        FIRSTLAST1(m.first,m.last) ); 
  }

  template <class T> 
  inline ConstSymBandMatrixView<T> HermBandMatrixViewOf(
      const GenBandMatrix<T>& m, UpLoType uplo, int nlo=-1)
  {
    if (nlo<0) nlo = (uplo==Upper?m.nhi():m.nlo());
    TMVAssert(m.colsize()==m.rowsize());
#ifdef XTEST_DEBUG
    TMVAssert(IsReal(T()) || m.diag().Imag().NormInf() == RealType(T)(0));
#endif
    return ConstSymBandMatrixView<T>(m.cptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),Herm,uplo,m.stor(),m.ct()); 
  }

  template <class T, IndexStyle I> 
  inline ConstSymBandMatrixView<T> HermBandMatrixViewOf(
      const ConstBandMatrixView<T,I>& m, UpLoType uplo, int nlo=-1)
  { 
    if (nlo<0) nlo = (uplo==Upper?m.nhi():m.nlo());
    TMVAssert(m.colsize()==m.rowsize());
#ifdef XTEST_DEBUG
    TMVAssert(IsReal(T()) || m.diag().Imag().NormInf() == RealType(T)(0));
#endif
    return ConstSymBandMatrixView<T,I>(m.cptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),Herm,uplo,m.stor(),m.ct()); 
  }

  template <class T, StorageType S, IndexStyle I> 
  inline ConstSymBandMatrixView<T,I> HermBandMatrixViewOf(
      const BandMatrix<T,S,I>& m, UpLoType uplo, int nlo=-1)
  {
    if (nlo<0) nlo = (uplo==Upper?m.nhi():m.nlo());
    TMVAssert(m.colsize()==m.rowsize());
#ifdef XTEST_DEBUG
    TMVAssert(IsReal(T()) || m.diag().Imag().NormInf() == RealType(T)(0));
#endif
    return ConstSymBandMatrixView<T,I>(m.cptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),Herm,uplo,m.stor(),m.ct()); 
  }

  template <class T, IndexStyle I> 
  inline SymBandMatrixView<T> HermBandMatrixViewOf(
      const BandMatrixView<T,I>& m, UpLoType uplo, int nlo=-1)
  { 
    if (nlo<0) nlo = (uplo==Upper?m.nhi():m.nlo());
    TMVAssert(m.colsize()==m.rowsize());
#ifdef XTEST_DEBUG
    TMVAssert(IsReal(T()) || m.diag().Imag().NormInf() == RealType(T)(0));
#endif
    return SymBandMatrixView<T,I>(m.ptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),Herm,uplo,m.stor(),m.ct()
        FIRSTLAST1(m.first,m.last) ); 
  }

  template <class T, StorageType S, IndexStyle I> 
  inline SymBandMatrixView<T,I> HermBandMatrixViewOf(
      BandMatrix<T,S,I>& m, UpLoType uplo, int nlo=-1)
  {
    if (nlo<0) nlo = (uplo==Upper?m.nhi():m.nlo());
    TMVAssert(m.colsize()==m.rowsize());
#ifdef XTEST_DEBUG
    TMVAssert(IsReal(T()) || m.diag().Imag().NormInf() == RealType(T)(0));
#endif
    return SymBandMatrixView<T,I>(m.ptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),Herm,uplo,m.stor(),m.ct()
        FIRSTLAST1(m.first,m.last) ); 
  }

  // From SymMatrix
  template <class T> 
  inline ConstSymBandMatrixView<T> SymBandMatrixViewOf(
      const GenSymMatrix<T>& m, int nlo)
  {
    return ConstSymBandMatrixView<T>(m.cptr(),m.size(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
        m.stor(),m.ct()); 
  }

  template <class T, IndexStyle I> 
  inline ConstSymBandMatrixView<T,I> SymBandMatrixViewOf(
      const ConstSymMatrixView<T,I>& m, int nlo)
  { 
    return ConstSymBandMatrixView<T,I>(m.cptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
        m.stor(),m.ct()); 
  }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline ConstSymBandMatrixView<T,I> SymBandMatrixViewOf(
      const SymMatrix<T,U,S,I>& m, int nlo)
  { 
    return ConstSymBandMatrixView<T,I>(m.cptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
        m.stor(),m.ct()); 
  }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline ConstSymBandMatrixView<T,I> SymBandMatrixViewOf(
      const HermMatrix<T,U,S,I>& m, int nlo)
  { 
    return ConstSymBandMatrixView<T,I>(m.cptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
        m.stor(),m.ct()); 
  }

  template <class T, IndexStyle I> 
  inline SymBandMatrixView<T,I> SymBandMatrixViewOf(
      const SymMatrixView<T,I>& m, int nlo)
  { 
    return SymBandMatrixView<T,I>(m.ptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
        m.stor(),m.ct() FIRSTLAST1(m.first,m.last) ); 
  }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline SymBandMatrixView<T,I> SymBandMatrixViewOf(
      SymMatrix<T,U,S,I>& m, int nlo)
  {
    return SymBandMatrixView<T,I>(m.ptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
        m.stor(),m.ct() FIRSTLAST1(m.first,m.last) ); 
  }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline SymBandMatrixView<T,I> SymBandMatrixViewOf(
      HermMatrix<T,U,S,I>& m, int nlo)
  {
    return SymBandMatrixView<T,I>(m.ptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
        m.stor(),m.ct() FIRSTLAST1(m.first,m.last) ); 
  }

  template <class T> 
  inline ConstSymBandMatrixView<T> HermBandMatrixViewOf(
      const GenSymMatrix<T>& m, int nlo)
  {
    TMVAssert(m.isherm());
    return ConstSymBandMatrixView<T>(m.cptr(),m.size(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
        m.stor(),m.ct()); 
  }

  template <class T, IndexStyle I> 
  inline ConstSymBandMatrixView<T,I> HermBandMatrixViewOf(
      const ConstSymMatrixView<T,I>& m, int nlo)
  { 
    TMVAssert(m.isherm());
    return ConstSymBandMatrixView<T,I>(m.cptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
        m.stor(),m.ct()); 
  }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline ConstSymBandMatrixView<T,I> HermBandMatrixViewOf(
      const SymMatrix<T,U,S,I>& m, int nlo)
  { 
    TMVAssert(m.isherm());
    return ConstSymBandMatrixView<T,I>(m.cptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
        m.stor(),m.ct()); 
  }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline ConstSymBandMatrixView<T,I> HermBandMatrixViewOf(
      const HermMatrix<T,U,S,I>& m, int nlo)
  { 
    TMVAssert(m.isherm());
    return ConstSymBandMatrixView<T,I>(m.cptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
        m.stor(),m.ct()); 
  }

  template <class T, IndexStyle I> 
  inline SymBandMatrixView<T,I> HermBandMatrixViewOf(
      const SymMatrixView<T,I>& m, int nlo)
  { 
    TMVAssert(m.isherm());
    return SymBandMatrixView<T,I>(m.ptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
        m.stor(),m.ct() FIRSTLAST1(m.first,m.last) ); 
  }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline SymBandMatrixView<T,I> HermBandMatrixViewOf(
      SymMatrix<T,U,S,I>& m, int nlo)
  {
    TMVAssert(m.isherm());
    return SymBandMatrixView<T,I>(m.ptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
        m.stor(),m.ct() FIRSTLAST1(m.first,m.last) ); 
  }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline SymBandMatrixView<T,I> HermBandMatrixViewOf(
      HermMatrix<T,U,S,I>& m, int nlo)
  {
    TMVAssert(m.isherm());
    return SymBandMatrixView<T,I>(m.ptr(),m.colsize(),nlo,
        m.stepi(),m.stepj(),m.stepi()+m.stepj(),m.sym(),m.uplo(),
        m.stor(),m.ct() FIRSTLAST1(m.first,m.last) ); 
  }

  // From ptr
  template <class T> 
  ConstSymBandMatrixView<T> SymBandMatrixViewOf(
      const T* vv, size_t size, int nlo, UpLoType uplo, StorageType stor);

  template <class T> 
  ConstSymBandMatrixView<T> HermBandMatrixViewOf(
      const T* vv, size_t size, int nlo, UpLoType uplo, StorageType stor);

  template <class T> 
  SymBandMatrixView<T> SymBandMatrixViewOf(
      T* vv, size_t size, int nlo, UpLoType uplo, StorageType stor);

  template <class T> 
  SymBandMatrixView<T> HermBandMatrixViewOf(
      T* vv, size_t size, int nlo, UpLoType uplo, StorageType stor);

  //
  // Swap Matrices
  //

  template <class T> 
  inline void Swap(
      const SymBandMatrixView<T>& m1, const SymBandMatrixView<T>& m2)
  { 
    TMVAssert(m1.size() == m2.size());
    TMVAssert(m1.nlo() == m2.nlo());
    TMVAssert(m1.issym() == m2.issym());
    TMVAssert(m1.isherm() == m2.isherm());
    Swap(m1.UpperBand(),m2.UpperBand()); 
  }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline void Swap(
      const SymBandMatrixView<T>& m1, SymBandMatrix<T,U,S,I>& m2)
  { Swap(m1,m2.View()); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline void Swap(
      SymBandMatrix<T,U,S,I>& m1, const SymBandMatrixView<T>& m2)
  { Swap(m1.View(),m2); }

  template <class T, UpLoType U1, StorageType S1, IndexStyle I1, UpLoType U2, StorageType S2, IndexStyle I2>
  inline void Swap(SymBandMatrix<T,U1,S1,I1>& m1, 
      SymBandMatrix<T,U2,S2,I2>& m2)
  { Swap(m1.View(),m2.View()); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline void Swap(
      const SymBandMatrixView<T>& m1, HermBandMatrix<T,U,S,I>& m2)
  { Swap(m1,m2.View()); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline void Swap(
      HermBandMatrix<T,U,S,I>& m1, const SymBandMatrixView<T>& m2)
  { Swap(m1.View(),m2); }

  template <class T, UpLoType U1, StorageType S1, IndexStyle I1, UpLoType U2, StorageType S2, IndexStyle I2>
  inline void Swap(HermBandMatrix<T,U1,S1,I1>& m1, 
      HermBandMatrix<T,U2,S2,I2>& m2)
  { Swap(m1.View(),m2.View()); }


  //
  // Views:
  //

  template <class T> 
  inline ConstSymBandMatrixView<T> Transpose(const GenSymBandMatrix<T>& m)
  { return m.Transpose(); }

  template <class T, IndexStyle I> 
  inline ConstSymBandMatrixView<T,I> Transpose(
      const ConstSymBandMatrixView<T,I>& m)
  { return m.Transpose(); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline ConstSymBandMatrixView<T,I> Transpose(
      const SymBandMatrix<T,U,S,I>& m)
  { return m.Transpose(); }

  template <class T, IndexStyle I> 
  inline SymBandMatrixView<T,I> Transpose(const SymBandMatrixView<T,I>& m)
  { return m.Transpose(); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline SymBandMatrixView<T,I> Transpose(SymBandMatrix<T,U,S,I>& m)
  { return m.Transpose(); }

  template <class T> 
  inline ConstSymBandMatrixView<T> Conjugate(const GenSymBandMatrix<T>& m)
  { return m.Conjugate(); }

  template <class T, IndexStyle I> 
  inline ConstSymBandMatrixView<T,I> Conjugate(
      const ConstSymBandMatrixView<T,I>& m)
  { return m.Conjugate(); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline ConstSymBandMatrixView<T,I> Conjugate(
      const SymBandMatrix<T,U,S,I>& m)
  { return m.Conjugate(); }

  template <class T, IndexStyle I> 
  inline SymBandMatrixView<T,I> Conjugate(const SymBandMatrixView<T,I>& m)
  { return m.Conjugate(); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline SymBandMatrixView<T,I> Conjugate(SymBandMatrix<T,U,S,I>& m)
  { return m.Conjugate(); }

  template <class T> 
  inline ConstSymBandMatrixView<T> Adjoint(const GenSymBandMatrix<T>& m)
  { return m.Adjoint(); }

  template <class T, IndexStyle I> 
  inline ConstSymBandMatrixView<T,I> Adjoint(
      const ConstSymBandMatrixView<T,I>& m)
  { return m.Adjoint(); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline ConstSymBandMatrixView<T,I> Adjoint(const SymBandMatrix<T,U,S,I>& m)
  { return m.Adjoint(); }

  template <class T, IndexStyle I> 
  inline SymBandMatrixView<T,I> Adjoint(const SymBandMatrixView<T,I>& m)
  { return m.Adjoint(); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline SymBandMatrixView<T,I> Adjoint(SymBandMatrix<T,U,S,I>& m)
  { return m.Adjoint(); }

  template <class T> 
  inline QuotXsB<T,T> Inverse(const GenSymBandMatrix<T>& m)
  { return m.Inverse(); }

  //
  // SymBandMatrix ==, != SymBandMatrix
  //

  template <class T1, class T2> 
  inline bool operator==(
      const GenSymBandMatrix<T1>& m1, const GenSymBandMatrix<T2>& m2)
  { return m1.UpperBand() == m2.UpperBand(); }
  template <class T1, class T2> 
  inline bool operator!=(
      const GenSymBandMatrix<T1>& m1, const GenSymBandMatrix<T2>& m2)
  { return !(m1 == m2); }


  //
  // I/O
  //

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  std::istream& operator>>(std::istream& is, 
      auto_ptr<SymBandMatrix<T,U,S,I> >& m);

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  std::istream& operator>>(std::istream& is,
      auto_ptr<HermBandMatrix<T,U,S,I> >& m);

  template <class T> 
  std::istream& operator>>(std::istream& is, const SymBandMatrixView<T>& m);

  template <class T, UpLoType U, StorageType S, IndexStyle I>
  inline std::istream& operator>>(std::istream& is, 
      SymBandMatrix<T,U,S,I>& m)
  { return is>>m.View(); }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  inline std::istream& operator>>(std::istream& is, 
      HermBandMatrix<T,U,S,I>& m)
  { return is>>m.View(); }

} // namespace tmv

#endif
