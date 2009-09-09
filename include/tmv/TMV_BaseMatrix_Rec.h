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

//-----------------------------------------------------------------------------
//
// This file defines the BaseMatrix_Rec and BaseMatrix_Rec_Mutable
// classes.  These are the base classes for the dense rectangular
// matrices, Matrix and SmallMatrix.
//
// BaseMatrix_Rec is the base class for all dense rectangular matrices.
// It adds the rest of the views described in TMV_Matrix.h like
// SubMatrix, SubVector, Rows, Cols, etc. that aren't defined for 
// a general BaseMatrix.
//
// BaseMatrix_Rec_Mutable is the base class for mutable dense 
// rectangular matrices and adds the non-const views along with the
// other modifying functions like TransposeSelf(), SetToIdentity(), etc.

#ifndef TMV_BaseMatrix_Rec_H
#define TMV_BaseMatrix_Rec_H

#include <sstream>
#include "TMV_BaseMatrix.h"

namespace tmv {

  // BaseMatrix_Rec is derived from BaseMatrix_Calc, and is used
  // for dense rectangular matrices.
  template <class M>
  class BaseMatrix_Rec;

  // BaseMatrix_Rec adds the following requirements to Traits<M>:
  //
  //  mstepi = the step size along column if known (else UNKNOWN)
  //  mstepj = the step size along row if known (else UNKNOWN)
  //  mdiagstep = the step size along row if known (else UNKNOWN)
  //
  //  const_row_type = return type from row(i) const
  //  const_col_type = return type from col(j) const
  //  const_diag_type = return type from diag() const
  //  const_row_range_type = return type from row(i,j1,j2) const
  //  const_col_range_type = return type from col(j,i1,i2) const
  //  const_diag_range_type = return type from diag(i) and diag(i,j1,j2) const
  //
  //  const_submatrix_type = return type from SubMatrix(i1,i2,j1,j2) const
  //  const_submatrix_step_type = return type from SubMatrix(i1,i2,j1,j2,istep,jstep) const
  //  const_subvector_type = return type from SubVector(i,j,istep,jstep,size) const
  //  const_colpair_type = return type from ColPair(j1,j2) const
  //  const_rowpair_type = return type from RowPair(i1,i2) const
  //  const_cols_type = return type from Cols(j1,j2) const
  //  const_rows_type = return type from Rows(i1,i2) const
  //
  //  const_uppertri_type = return type from UpperTri() const
  //  const_unit_uppertri_type = return type from UnitUpperTri() const
  //  const_lowertri_type = return type from LowerTri() const
  //  const_unit_lowertri_type = return type from UnitLowerTri() const
  //  const_linearview_type = return type from LinearView() const
  //  


  // BaseMatrix_Rec_Mutable is derived from both BaseMatrix_Rec and
  // BaseMatrix_Mutable, and is used for mutable, dense, rectangular matrices.
  template <class M>
  class BaseMatrix_Rec_Mutable;

  // BaseMatrix_Rec_Mutable adds:
  //
  //  col_type = the return type of col(i) 
  //  row_type = the return type of row(j) 
  //  diag_type = return type from diag() 
  //  row_range_type = return type from row(i,j1,j2) 
  //  col_range_type = return type from col(j,i1,i2) 
  //  diag_range_type = return type from diag(i) and diag(i,j1,j2) 
  //
  //  submatrix_type = return type from SubMatrix(i1,i2,j1,j2) 
  //  submatrix_step_type = return type from SubMatrix(i1,i2,j1,j2,istep,jstep) 
  //  subvector_type = return type from SubVector(i,j,istep,jstep,size)
  //  colpair_type = return type from ColPair(j1,j2) 
  //  rowpair_type = return type from RowPair(i1,i2) 
  //  cols_type = return type from Cols(j1,j2) 
  //  rows_type = return type from Rows(i1,i2) 
  //
  //  uppertri_type = return type from UpperTri() 
  //  unit_uppertri_type = return type from UnitUpperTri() 
  //  lowertri_type = return type from LowerTri() 
  //  unit_lowertri_type = return type from UnitLowerTri() 
  //  realview_type = return type from Real()
  //  imagview_type = return type from Imag() 
  //  nonconj_type = return type from NonConj() 
  //  linearview_type = return type from LinearView() 
  //  

  // The following all derive from BaseMatrix_Rec or BaseMatrix_Rec_Mutable.
  // See TMV_Matrix.h and TMV_SmallMatrix.h for their definitions:
  template <class T, StorageType S=ColMajor, IndexStyle I=CStyle>
  class Matrix;
  template <class T, int Si=UNKNOWN, int Sj=UNKNOWN, bool C=false, IndexStyle I=CStyle>
  class ConstMatrixView;
  template <class T, int Si=UNKNOWN, int Sj=UNKNOWN, bool C=false, IndexStyle I=CStyle>
  class MatrixView;
  template <class T, int M, int N, StorageType S=ColMajor, IndexStyle I=CStyle>
  class SmallMatrix;
  template <class T, int M, int N, int Si, int Sj, bool C=false, IndexStyle I=CStyle>
  class ConstSmallMatrixView;
  template <class T, int M, int N, int Si, int Sj, bool C=false, IndexStyle I=CStyle>
  class SmallMatrixView;

  // These are effectively aliases for I = FortranStyle
  // so you don' thave to write all the Si,Sj,C values if you are using
  // the defaults when you want to spectify FortranStyle.
  template <class T, StorageType S=ColMajor>
  class MatrixF;
  template <class T, int Si=UNKNOWN, int Sj=UNKNOWN, bool C=false>
  class ConstMatrixViewF;
  template <class T, int Si=UNKNOWN, int Sj=UNKNOWN, bool C=false>
  class MatrixViewF;
  template <class T, int M, int N, StorageType S=ColMajor>
  class SmallMatrixF;
  template <class T, int M, int N, int Si, int Sj, bool C=false>
  class ConstSmallMatrixViewF;
  template <class T, int M, int N, int Si, int Sj, bool C=false>
  class SmallMatrixViewF;


  //
  // Helper functions and values:
  //

  // Specify ExactSameStorage for rectangular matrices:
  template <class M>
  inline bool ExactSameStorage(
      const BaseMatrix_Rec<M>& m1, const BaseMatrix_Rec<M>& m2)
  { return m1.stepi() == m2.stepi() && m1.stepj() == m2.stepj(); }

  // These helper functions check the validity of indices according
  // to whether the matrix uses CStyle or FortranStyle indexing.
  // They also update the indices to be consistent with CStyle.
  template <bool mfort>
  inline void CheckDiagIndex(int& i, int m, int n)
  { // CStyle or FortranStyle
    TMVAssert(i >= -m && "negative diag index must be <= nrows");
    TMVAssert(i <= n && "positive diag index must be <= ncols");
  }
  template <bool mfort>
  inline void CheckDiagIndex(int& i, int& j1, int& j2, int m, int n)
  { // CStyle
    TMVAssert(i >= -m && "negative diag index must be <= nrows");
    TMVAssert(i <= n && "positive diag index must be <= ncols");
    TMVAssert(j1 >= 0 && "first element must be in matrix");
    TMVAssert(j2 >= j1 && "range must have a non-negative number of elements");
    // if i >= 0, i + j2 <= n  and j2 <= m
    TMVAssert((i < 0 || i+j2 <= n) && "last element must be in matrix");
    TMVAssert((i < 0 || j2 <= m) && "last element must be in matrix");
    // if i <= 0, -i + j2 <= m  and j2 <= n
    TMVAssert((i > 0 || -i+j2 <= m) && "last element must be in matrix");
    TMVAssert((i > 0 || j2 <= n) && "last element must be in matrix");
  }
  template <bool mfort>
  inline void CheckRowRange(int& i1, int i2, int m)
  { // CStyle
    TMVAssert(i1 >= 0 && "first row must be in matrix");
    TMVAssert(i2 <= m && "last row must be in matrix");
    TMVAssert(i2 >= i1 && "range must have a non-negative number of rows");
  }
  template <bool mfort>
  inline void CheckColRange(int& j1, int j2, int n)
  { // CStyle
    TMVAssert(j1 >= 0 && "first column must be in matrix");
    TMVAssert(j2 <= n && "last column must be in matrix");
    TMVAssert(j2 >= j1 && "range must have a non-negative number of columns");
  }
  template <bool mfort>
  inline void CheckRowRange(int& i1, int& i2, int istep, int m)
  { // CStyle
    TMVAssert(istep != 0 && "istep cannot be 0");
    TMVAssert(((i1 >= 0 && i1 < m) || i1==i2) && 
        "first row must be in matrix");
    TMVAssert(((i2-istep >= 0 && i2-istep < m) || i1==i2) && 
        "last row must be in matrix");
    TMVAssert((i2-i1) % istep == 0 &&
        "row range must be an integral multiple of istep");
    TMVAssert((i2-i1) / istep >= 0 &&
        "must have a non-negative number of rows");
  }
  template <bool mfort>
  inline void CheckColRange(int& j1, int& j2, int jstep, int n)
  { // CStyle
    TMVAssert(jstep != 0 && "jstep cannot be 0");
    TMVAssert(((j1 >= 0 && j1 < n) || j1==j2) && 
        "first column must be in matrix");
    TMVAssert(((j2-jstep >= 0 && j2-jstep < n) || j1==j2) &&
        "last column must be in matrix");
    TMVAssert((j2-j1) % jstep == 0 &&
        "column range must be an integral multiple of jstep");
    TMVAssert((j2-j1) / jstep >= 0 &&
        "must have a non-negative number of columns");
  }
  template <bool mfort>
  inline void CheckMatSubVector(int& i, int& j, int istep, int jstep, 
      int size, int m, int n)
  { // CStyle
    TMVAssert(!(istep == 0 && jstep == 0) && 
        "istep and jstep cannot both be 0");
    TMVAssert(((i >= 0 && i < m) || size==0) && 
        "first element must be in matrix");
    TMVAssert(((j >= 0 && j < n) || size==0) && 
        "first element must be in matrix");
    TMVAssert(((i+istep*(size-1) >= 0 && i+istep*(size-1) < m) || size==0) && 
        "last element must be in matrix");
    TMVAssert(((j+jstep*(size-1) >= 0 && j+jstep*(size-1) < n) || size==0) && 
        "last element must be in matrix");
  }
  template <>
  inline void CheckDiagIndex<true>(int& i, int& j1, int& j2, int m, int n)
  { // FortranStyle
    TMVAssert(i >= -m && "negative diag index must be <= nrows");
    TMVAssert(i <= n && "positive diag index must be <= ncols");
    TMVAssert(j1 >= 1 && "first element must be in matrix");
    TMVAssert(j2 >= j1 && "range must have a positive number of elements");
    TMVAssert((i < 0 || i+j2 <= n) && "last element must be in matrix");
    TMVAssert((i < 0 || j2 <= m) && "last element must be in matrix");
    TMVAssert((i > 0 || -i+j2 <= m) && "last element must be in matrix");
    TMVAssert((i > 0 || j2 <= n) && "last element must be in matrix");
    --j1;
  }
  template <>
  inline void CheckRowRange<true>(int& i1, int i2, int m)
  { // FortranStyle
    TMVAssert(i1 >= 1 && "first row must be in matrix");
    TMVAssert(i2 <= m && "last row must be in matrix");
    TMVAssert(i2 >= i1 && "range must have a positive number of rows");
    --i1;
  }
  template <>
  inline void CheckColRange<true>(int& j1, int j2, int n)
  { // FortranStyle
    TMVAssert(j1 >= 1 && "first column must be in matrix");
    TMVAssert(j2 <= n && "last column must be in matrix");
    TMVAssert(j2 >= j1 && "range must have a positive number of columns");
    --j1;
  }
  template <>
  inline void CheckRowRange<true>(int& i1, int& i2, int istep, int m)
  { // FortranStyle
    TMVAssert(istep != 0 && "istep cannot be 0");
    TMVAssert(i1 >= 1 && i1 <= m && "first row must be in matrix");
    TMVAssert(i2 >= 1 && i2 <= m && "last row must be in matrix");
    TMVAssert((i2-i1) % istep == 0 &&
        "row range must be an integral multiple of istep");
    TMVAssert((i2-i1) / istep >= 0 &&
        "must have a positive number of rows");
    --i1; i2 += istep-1;
  }
  template <>
  inline void CheckColRange<true>(int& j1, int& j2, int jstep, int n)
  { // FortranStyle
    TMVAssert(jstep != 0 && "jstep cannot be 0");
    TMVAssert(j1 >= 1 && j1 <= n && "first column must be in matrix");
    TMVAssert(j2 >= 1 && j2 <= n && "last column must be in matrix");
    TMVAssert((j2-j1) % jstep == 0 &&
        "column range must be an integral multiple of jstep");
    TMVAssert((j2-j1) / jstep >= 0 &&
        "must have a positive number of columns");
    --j1; j2 += jstep-1;
  }
  template <>
  inline void CheckMatSubVector<true>(int& i, int& j, int istep, int jstep, 
      int size, int m, int n)
  { // FortranStyle
    TMVAssert(!(istep == 0 && jstep == 0) && 
        "istep and jstep cannot both be 0");
    TMVAssert(i >= 1 && i <= m && "first element must be in matrix");
    TMVAssert(j >= 1 && j <= n && "first element must be in matrix");
    TMVAssert(i+istep*(size-1) >= 1 && i+istep*(size-1) <= m && 
        "last element must be in matrix");
    TMVAssert(j+jstep*(size-1) >= 1 && j+jstep*(size-1) <= n && 
        "last element must be in matrix");
    --i; --j;
  }

  template <class T, int cs, int rs, bool rm, bool fort>
  struct MCopyHelper<T,Rec,cs,rs,rm,fort>
  { typedef SmallMatrix<T,cs,rs,rm?RowMajor:ColMajor,fort?FortranStyle:CStyle> type; };
  template <class T, int rs, bool rm, bool fort>
  struct MCopyHelper<T,Rec,UNKNOWN,rs,rm,fort>
  { typedef Matrix<T,rm?RowMajor:ColMajor,fort?FortranStyle:CStyle> type; };
  template <class T, int cs, bool rm, bool fort>
  struct MCopyHelper<T,Rec,cs,UNKNOWN,rm,fort>
  { typedef Matrix<T,rm?RowMajor:ColMajor,fort?FortranStyle:CStyle> type; };
  template <class T, bool rm, bool fort>
  struct MCopyHelper<T,Rec,UNKNOWN,UNKNOWN,rm,fort>
  { typedef Matrix<T,rm?RowMajor:ColMajor,fort?FortranStyle:CStyle> type; };

  // A quick auxilliary function for CanLinearize.
  // (It only accesses the steps that are unknown at compile time.)
  template <bool canlin, int Si, int Sj, class M>
  struct AuxCanLinearize;

  template <int Si, int Sj, class M>
  struct AuxCanLinearize<true,Si,Sj,M>
  {
    static inline bool ok(const M& m) 
    { return true; }
  };


  template <int Si, int Sj, class M>
  struct AuxCanLinearize<false,Si,Sj,M>
  {
    static inline bool ok(const M& m) 
    {
      return (
          (m.stepi() == 1 && m.stepj() == m.colsize()) ||
          (m.stepj() == 1 && m.stepi() == m.rowsize()) );
    }
  };
  template <int Sj, class M>
  struct AuxCanLinearize<false,1,Sj,M>
  { static inline bool ok(const M& m) { return m.stepj() == m.colsize(); } };
  template <int Si, class M>
  struct AuxCanLinearize<false,Si,1,M>
  { static inline bool ok(const M& m) { return m.stepi() == m.rowsize(); } };
  template <class M>
  struct AuxCanLinearize<false,1,1,M>
  { static inline bool ok(const M& m) { return true; } };

  // Defined in TMV_CopyM.h
  template <class M1, class M2>
  inline void Copy(
      const BaseMatrix_Rec<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2);

  // Defined in TMV_SwapM.h
  template <class M1, class M2>
  inline void Swap(
      BaseMatrix_Rec_Mutable<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2);
  template <class M>
  inline void TransposeSelf(BaseMatrix_Rec_Mutable<M>& m);

  // Defined in TMV_Permute.h
  template <class M>
  inline void PermuteRows(BaseMatrix_Rec_Mutable<M>& m,
      const int*const p, int i1, int i2);
  template <class M>
  inline void ReversePermuteRows(BaseMatrix_Rec_Mutable<M>& m,
      const int*const p, int i1, int i2);

  // Defined in TMV_NormM.h
  template <class M>
  inline typename M::real_type NormSq(const BaseMatrix_Rec<M>& m);
  template <class M>
  inline typename M::real_type NormSq(const BaseMatrix_Rec<M>& m, 
      const typename M::real_type scale);
  template <class M>
  inline typename M::real_type NormF(const BaseMatrix_Rec<M>& m);
  template <class M>
  inline typename M::real_type Norm1(const BaseMatrix_Rec<M>& m);
  template <class M>
  inline typename M::real_type NormInf(const BaseMatrix_Rec<M>& m);
  template <class M>
  inline typename M::value_type SumElements(const BaseMatrix_Rec<M>& m);
  template <class M>
  inline typename M::real_type SumAbsElements(const BaseMatrix_Rec<M>& m);
  template <class M>
  inline typename M::real_type MaxAbsElement(const BaseMatrix_Rec<M>& m);

  // Defined in TMV_MatrixIO.h
  template <class M>
  inline void Read(std::istream& is, BaseMatrix_Rec_Mutable<M>& m);

  // Defined below:
  template <class M>
  inline void Zero(BaseMatrix_Rec_Mutable<M>& m);
  template <class M, class RT>
  inline void Clip(BaseMatrix_Rec_Mutable<M>& m, const RT& thresh);
  template <class M, class T>
  inline void SetAllTo(BaseMatrix_Rec_Mutable<M>& m, const T& val);
  template <class M, class T>
  inline void AddToAll(BaseMatrix_Rec_Mutable<M>& m, const T& val);
  template <class M>
  inline void ConjugateSelf(BaseMatrix_Rec_Mutable<M>& m);
  template <class M, class F>
  inline void ApplyToAll(BaseMatrix_Rec_Mutable<M>& m, const F& f);


#if 0
  // Defined in TMV_MatrixDiv.h
  template <class M>
  inline typename M::value_type Det(const BaseMatrix_Rec<M>& m);
  template <class M>
  inline typename M::real_type LogDet(const BaseMatrix_Rec<M>& m, 
      typename M::value_type* sign=0);
  template <class M, class M2>
  inline void DoInverse(const BaseMatrix_Rec<M>& m, BaseMatrix_Rec_Mutable<M2>& m2);
  template <class M, class M2>
  inline void DoInverseATA(const BaseMatrix_Rec<M>& m, BaseMatrix_Rec_Mutable<M2>& m2);
  template <class M>
  inline typename M::real_type Norm2(const BaseMatrix_Rec<M>& m);
  template <class M>
  inline typename M::real_type Condition(const BaseMatrix_Rec<M>& m);
#endif

  template <class M> 
  class BaseMatrix_Rec : 
    public BaseMatrix_Calc<M>
  {
  public:

    typedef M type;

    typedef typename Traits<type>::value_type value_type;
    typedef typename Traits<type>::calc_type calc_type;
    typedef typename Traits<type>::eval_type eval_type;
    typedef typename Traits<type>::copy_type copy_type;

    enum { mcolsize = Traits<type>::mcolsize };
    enum { mrowsize = Traits<type>::mrowsize };
    enum { mshape = Traits<type>::mshape };
    enum { mfort = Traits<type>::mfort };
    enum { mrowmajor = Traits<type>::mrowmajor }; 
    enum { mcolmajor = Traits<type>::mcolmajor }; 
    enum { mstor = Traits<type>::mstor };
    enum { mcalc = Traits<type>::mcalc };
    enum { mstepi = Traits<type>::mstepi };
    enum { mstepj = Traits<type>::mstepj };
    enum { mdiagstep = Traits<type>::mdiagstep };
    enum { mconj = Traits<type>::mconj };
    enum { mcanlin = Traits<type>::mcanlin };

    typedef typename Traits<value_type>::real_type real_type;
    typedef typename Traits<value_type>::complex_type complex_type;
    enum { misreal = Traits<value_type>::isreal };
    enum { miscomplex = Traits<value_type>::iscomplex };

    typedef typename Traits<type>::const_row_type const_row_type;
    typedef typename Traits<type>::const_row_range_type const_row_range_type;
    typedef typename Traits<type>::const_col_type const_col_type;
    typedef typename Traits<type>::const_col_range_type const_col_range_type;
    typedef typename Traits<type>::const_diag_type const_diag_type;
    typedef typename Traits<type>::const_diag_range_type const_diag_range_type;

    typedef typename Traits<type>::const_submatrix_type const_submatrix_type;
    typedef typename Traits<type>::const_submatrix_step_type const_submatrix_step_type;
    typedef typename Traits<type>::const_subvector_type const_subvector_type;
    typedef typename Traits<type>::const_colpair_type const_colpair_type;
    typedef typename Traits<type>::const_rowpair_type const_rowpair_type;
    typedef typename Traits<type>::const_cols_type const_cols_type;
    typedef typename Traits<type>::const_rows_type const_rows_type;

    typedef typename Traits<type>::const_view_type const_view_type;
    typedef typename Traits<type>::const_cview_type const_cview_type;
    typedef typename Traits<type>::const_fview_type const_fview_type;
    typedef typename Traits<type>::const_xview_type const_xview_type;
    typedef typename Traits<type>::const_cmview_type const_cmview_type;
    typedef typename Traits<type>::const_rmview_type const_rmview_type;
    typedef typename Traits<type>::const_transpose_type const_transpose_type;
    typedef typename Traits<type>::const_conjugate_type const_conjugate_type;
    typedef typename Traits<type>::const_adjoint_type const_adjoint_type;
    typedef typename Traits<type>::const_uppertri_type const_uppertri_type;
    typedef typename Traits<type>::const_unit_uppertri_type const_unit_uppertri_type;
    typedef typename Traits<type>::const_lowertri_type const_lowertri_type;
    typedef typename Traits<type>::const_unit_lowertri_type const_unit_lowertri_type;
    typedef typename Traits<type>::const_linearview_type const_linearview_type;
    typedef typename Traits<type>::const_realview_type const_realview_type;
    typedef typename Traits<type>::const_imagview_type const_imagview_type;
    typedef typename Traits<type>::const_nonconj_type const_nonconj_type;
    typedef typename Traits<type>::nonconst_type nonconst_type;



    //
    // Constructor
    //

    inline BaseMatrix_Rec() {}
    inline BaseMatrix_Rec(const BaseMatrix_Rec<M>&) {}
    inline ~BaseMatrix_Rec() {}

  private:
    void operator=(const BaseMatrix_Rec<M>&);
  public:


    //
    // Access 
    //

    // The get_ routines always use CStyle indexing.
    inline const_row_type get_row(int i) const
    { return const_row_type(cptr()+i*stepi(),rowsize(),stepj()); }

    inline const_col_type get_col(int j) const
    { return const_col_type(cptr()+j*stepj(),colsize(),stepi()); }

    inline const_row_range_type get_row(int i, int j1, int j2) const
    { return const_row_range_type(cptr()+i*stepi()+j1*stepj(),j2-j1,stepj()); }

    inline const_col_range_type get_col(int j, int i1, int i2) const
    { return const_col_range_type(cptr()+j*stepj()+i1*stepi(),i2-i1,stepi()); }

    inline const_diag_range_type get_diag(int i) const
    {
      return const_diag_range_type(
          cptr() + (i<0?(-i*stepi()):(i*stepj())),
          (i<0? TMV_MIN(colsize()+i,rowsize()) : TMV_MIN(colsize(),rowsize()-i)),
          diagstep());
    }

    inline const_diag_range_type get_diag(int i, int j1, int j2) const
    {
      return const_diag_range_type(
          cptr() + (i<0?(-i*stepi()):(i*stepj())) + j1*diagstep(),
          j2-j1, diagstep());
    }


    // The regular versions respect the indexing style for i and j:
    inline const_row_type row(int i) const
    {
      CheckRowIndex<mfort>(i,colsize());
      return get_row(i);
    }

    inline const_row_type operator[](int i) const
    { return row(i); }

    inline const_col_type col(int j) const
    {
      CheckColIndex<mfort>(j,rowsize());
      return get_col(j);
    }

    inline const_row_range_type row(int i, int j1, int j2) const
    {
      CheckRowIndex<mfort>(i,colsize());
      CheckColRange<mfort>(j1,j2,rowsize());
      return get_row(i,j1,j2); 
    }

    inline const_col_range_type col(int j, int i1, int i2) const
    {
      CheckColIndex<mfort>(j,rowsize());
      CheckRowRange<mfort>(i1,i2,colsize());
      return get_col(j,i1,i2); 
    }

    // No need for a get_ routine for diag()
    inline const_diag_type diag() const
    { return const_diag_type(cptr(),TMV_MIN(colsize(),rowsize()),diagstep()); }

    inline const_diag_range_type diag(int i) const
    {
      CheckDiagIndex<mfort>(i,colsize(),rowsize());
      return get_diag(i);
    }

    inline const_diag_range_type diag(int i, int j1, int j2) const
    {
      CheckDiagIndex<mfort>(i,j1,j2,colsize(),rowsize());
      return get_diag(i,j1,j2);
    }


    //
    // Functions
    //

    inline value_type SumElements() const
    { return tmv::SumElements(CView()); }

    inline real_type SumAbsElements() const
    { return tmv::SumAbsElements(CView()); }

    inline real_type MaxAbsElement() const
    { return tmv::MaxAbsElement(CView()); }

    inline real_type NormSq() const
    { return tmv::NormSq(CView()); }

    inline real_type NormSq(const real_type scale) const
    { return tmv::NormSq(CView(),scale); }

    inline real_type NormF() const 
    { return tmv::NormF(CView()); }

    inline real_type Norm() const
    { return NormF(); }

    inline real_type Norm1() const
    { return tmv::Norm1(CView()); }

    inline real_type NormInf() const
    { return tmv::NormInf(CView()); }


    
    //
    // SubMatrix, etc.
    //

    // These versions always uses CStyle
    inline const_submatrix_type CSubMatrix(
        int i1, int i2, int j1, int j2) const
    {
      return const_submatrix_type(
          cptr()+i1*stepi()+j1*stepj(), i2-i1, j2-j1, stepi(), stepj()); 
    }

    inline const_submatrix_step_type CSubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) const
    {
      return const_submatrix_step_type(
          cptr()+i1*stepi()+j1*stepj(), (i2-i1)/istep, (j2-j1)/jstep,
          istep*stepi(), jstep*stepj());
    }

    inline const_subvector_type CSubVector(
        int i, int j, int istep, int jstep, int s) const
    {
      return const_subvector_type(
          cptr()+i*stepi()+j*stepj(), s, istep*stepi() + jstep*stepj());
    }

    inline const_colpair_type CColPair(int j1, int j2) const
    {
      return const_colpair_type(
          cptr()+j1*stepj(), colsize(), 2, stepi(), (j2-j1)*stepj());
    }

    inline const_rowpair_type CRowPair(int i1, int i2) const
    {
      return const_rowpair_type(
          cptr()+i1*stepi(), 2, rowsize(), (i2-i1)*stepi(), stepj());
    }

    inline const_cols_type CCols(int j1, int j2) const
    {
      return const_cols_type(
          cptr()+j1*stepj(), colsize(), j2-j1, stepi(), stepj());
    }

    inline const_rows_type CRows(int i1, int i2) const
    {
      return const_rows_type(
          cptr()+i1*stepi(), i2-i1, rowsize(), stepi(), stepj());
    }


    // These check the indices according the the indexing style being
    // used, and then calls the above CStyle versions.
    inline const_submatrix_type SubMatrix(
        int i1, int i2, int j1, int j2) const
    {
      CheckRowRange<mfort>(i1,i2,colsize());
      CheckColRange<mfort>(j1,j2,rowsize());
      return CSubMatrix(i1,i2,j1,j2);
    }

    inline const_submatrix_step_type SubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) const
    {
      CheckRowRange<mfort>(i1,i2,istep,colsize());
      CheckColRange<mfort>(j1,j2,jstep,rowsize());
      return CSubMatrix(i1,i2,j1,j2,istep,jstep);
    }

    inline const_subvector_type SubVector(
        int i, int j, int istep, int jstep, int s) const
    {
      CheckMatSubVector<mfort>(i,j,istep,jstep,s,colsize(),rowsize());
      return CSubVector(i,j,istep,jstep,s); 
    }

    inline const_colpair_type ColPair(int j1, int j2) const
    {
      CheckColIndex<mfort>(j1,rowsize());
      CheckColIndex<mfort>(j2,rowsize());
      return CColPair(j1,j2);
    }

    inline const_rowpair_type RowPair(int i1, int i2) const
    {
      CheckRowIndex<mfort>(i1,colsize());
      CheckRowIndex<mfort>(i2,colsize());
      return CRowPair(i1,i2);
    }

    inline const_cols_type Cols(int j1, int j2) const
    {
      CheckColRange<mfort>(j1,j2,rowsize());
      return CCols(j1,j2);
    }

    inline const_rows_type Rows(int i1, int i2) const
    {
      CheckRowRange<mfort>(i1,i2,colsize());
      return CRows(i1,i2);
    }


    //
    // Views
    //

    inline const_view_type View() const
    { return const_view_type(cptr(),colsize(),rowsize(),stepi(),stepj()); }

    inline const_cview_type CView() const
    { return View(); }

    inline const_fview_type FView() const
    { return View(); }

    inline const_xview_type XView() const
    { return View(); }

    inline const_cmview_type CMView() const
    { return View(); }

    inline const_rmview_type RMView() const
    { return View(); }

    inline const_transpose_type Transpose() const
    { return const_transpose_type(cptr(),rowsize(),colsize(),stepj(),stepi()); }

    inline const_conjugate_type Conjugate() const
    { return const_conjugate_type(cptr(),colsize(),rowsize(),stepi(),stepj()); }

    inline const_adjoint_type Adjoint() const
    { return const_adjoint_type(cptr(),rowsize(),colsize(),stepj(),stepi()); }

    inline const_uppertri_type UpperTri() const
    { return const_uppertri_type(cptr(),rowsize(),stepi(),stepj()); }

    inline const_unit_uppertri_type UnitUpperTri() const
    { return const_unit_uppertri_type(cptr(),rowsize(),stepi(),stepj()); }

    inline const_lowertri_type LowerTri() const
    { return const_lowertri_type(cptr(),colsize(),stepi(),stepj()); }

    inline const_unit_lowertri_type UnitLowerTri() const
    { return const_unit_lowertri_type(cptr(),colsize(),stepi(),stepj()); }

    inline const_linearview_type LinearView() const
    {
      TMVAssert(
          (stepi() == 1 && stepj() == colsize()) ||
          (stepj() == 1 && stepi() == rowsize()) );
      return const_linearview_type(cptr(),ls(),1);
    }

    inline const_realview_type Real() const
    {
      return const_realview_type(
          reinterpret_cast<const real_type*>(cptr()), colsize(), rowsize(),
          misreal ? stepi() : 2*stepi(), misreal ? stepj() : 2*stepj());
    }

    inline const_imagview_type Imag() const
    {
      TMVStaticAssert(miscomplex);
      return const_imagview_type(
          reinterpret_cast<const real_type*>(cptr())+1, colsize(), rowsize(),
          misreal ? stepi() : 2*stepi(), misreal ? stepj() : 2*stepj());
    }

    inline const_nonconj_type NonConj() const
    { return const_nonconj_type(cptr(),colsize(),rowsize(),stepi(),stepj()); }

    inline nonconst_type NonConst() const
    {
      return nonconst_type(const_cast<value_type*>(cptr()),
          colsize(),rowsize(),stepi(),stepj());
    }



    //
    // Auxilliary routines
    //

    template <class M2>
    inline void AssignTo(BaseMatrix_Mutable<M2>& m2) const
    {
      TMVStaticAssert((ShapeTraits2<mshape,M2::mshape>::assignable));
      tmv::Copy(mat(),m2.mat()); 
    }

    inline const type& mat() const
    { return *static_cast<const type*>(this); }

    inline int diagstep() const 
    { return mdiagstep == UNKNOWN ? stepi() + stepj() : mdiagstep; }
    inline bool isconj() const { return mconj; }
    inline bool CanLinearize() const 
    { return (AuxCanLinearize<mcanlin,mcolmajor,mrowmajor,M>::ok(mat())); }

    // Note that these last functions need to be defined in a more derived
    // class than this, or an infinite loop will result when compiling.
    // Also, cref from BaseMatrix.

    inline size_t colsize() const { return mat().colsize(); }
    inline size_t rowsize() const { return mat().rowsize(); }
    inline size_t ls() const { return mat().ls(); }
    inline int stepi() const { return mat().stepi(); }
    inline int stepj() const { return mat().stepj(); }
    inline bool isrm() const { return mat().isrm(); }
    inline bool iscm() const { return mat().iscm(); }

    inline const value_type* cptr() const { return mat().cptr(); }

  }; // BaseMatrix_Rec

  template <class M> 
  class BaseMatrix_Rec_Mutable : 
    public BaseMatrix_Rec<M>,
    public BaseMatrix_Mutable<M>
  {
  public:

    typedef M type;
    typedef BaseMatrix_Rec<M> base_rec;

    typedef typename Traits<type>::value_type value_type;
    typedef typename Traits<type>::calc_type calc_type;
    typedef typename Traits<type>::eval_type eval_type;
    typedef typename Traits<type>::copy_type copy_type;

    enum { mcolsize = Traits<type>::mcolsize };
    enum { mrowsize = Traits<type>::mrowsize };
    enum { mshape = Traits<type>::mshape };
    enum { mfort = Traits<type>::mfort };
    enum { mcalc = Traits<type>::mcalc };
    enum { mrowmajor = Traits<type>::mrowmajor }; 
    enum { mcolmajor = Traits<type>::mcolmajor }; 
    enum { mstor = Traits<type>::mstor };
    enum { mstepi = Traits<type>::mstepi };
    enum { mstepj = Traits<type>::mstepj };
    enum { mdiagstep = Traits<type>::mdiagstep };
    enum { mconj = Traits<type>::mconj };
    enum { mcanlin = Traits<type>::mcanlin };

    typedef typename Traits<value_type>::real_type real_type;
    typedef typename Traits<value_type>::complex_type complex_type;
    enum { misreal = Traits<value_type>::isreal };
    enum { miscomplex = Traits<value_type>::iscomplex };

    typedef typename Traits<type>::const_row_type const_row_type;
    typedef typename Traits<type>::const_row_range_type const_row_range_type;
    typedef typename Traits<type>::const_col_type const_col_type;
    typedef typename Traits<type>::const_col_range_type const_col_range_type;
    typedef typename Traits<type>::const_diag_type const_diag_type;
    typedef typename Traits<type>::const_diag_range_type const_diag_range_type;

    typedef typename Traits<type>::const_submatrix_type const_submatrix_type;
    typedef typename Traits<type>::const_submatrix_step_type const_submatrix_step_type;
    typedef typename Traits<type>::const_subvector_type const_subvector_type;
    typedef typename Traits<type>::const_colpair_type const_colpair_type;
    typedef typename Traits<type>::const_rowpair_type const_rowpair_type;
    typedef typename Traits<type>::const_cols_type const_cols_type;
    typedef typename Traits<type>::const_rows_type const_rows_type;

    typedef typename Traits<type>::const_view_type const_view_type;
    typedef typename Traits<type>::const_cview_type const_cview_type;
    typedef typename Traits<type>::const_fview_type const_fview_type;
    typedef typename Traits<type>::const_xview_type const_xview_type;
    typedef typename Traits<type>::const_cmview_type const_cmview_type;
    typedef typename Traits<type>::const_rmview_type const_rmview_type;
    typedef typename Traits<type>::const_transpose_type const_transpose_type;
    typedef typename Traits<type>::const_conjugate_type const_conjugate_type;
    typedef typename Traits<type>::const_adjoint_type const_adjoint_type;
    typedef typename Traits<type>::const_uppertri_type const_uppertri_type;
    typedef typename Traits<type>::const_unit_uppertri_type const_unit_uppertri_type;
    typedef typename Traits<type>::const_lowertri_type const_lowertri_type;
    typedef typename Traits<type>::const_unit_lowertri_type const_unit_lowertri_type;
    typedef typename Traits<type>::const_linearview_type const_linearview_type;
    typedef typename Traits<type>::const_realview_type const_realview_type;
    typedef typename Traits<type>::const_imagview_type const_imagview_type;
    typedef typename Traits<type>::const_nonconj_type const_nonconj_type;

    typedef typename Traits<type>::row_type row_type;
    typedef typename Traits<type>::row_range_type row_range_type;
    typedef typename Traits<type>::col_type col_type;
    typedef typename Traits<type>::col_range_type col_range_type;
    typedef typename Traits<type>::diag_type diag_type;
    typedef typename Traits<type>::diag_range_type diag_range_type;

    typedef typename Traits<type>::submatrix_type submatrix_type;
    typedef typename Traits<type>::submatrix_step_type submatrix_step_type;
    typedef typename Traits<type>::subvector_type subvector_type;
    typedef typename Traits<type>::colpair_type colpair_type;
    typedef typename Traits<type>::rowpair_type rowpair_type;
    typedef typename Traits<type>::cols_type cols_type;
    typedef typename Traits<type>::rows_type rows_type;

    typedef typename Traits<type>::view_type view_type;
    typedef typename Traits<type>::cview_type cview_type;
    typedef typename Traits<type>::fview_type fview_type;
    typedef typename Traits<type>::xview_type xview_type;
    typedef typename Traits<type>::cmview_type cmview_type;
    typedef typename Traits<type>::rmview_type rmview_type;
    typedef typename Traits<type>::transpose_type transpose_type;
    typedef typename Traits<type>::conjugate_type conjugate_type;
    typedef typename Traits<type>::adjoint_type adjoint_type;
    typedef typename Traits<type>::uppertri_type uppertri_type;
    typedef typename Traits<type>::unit_uppertri_type unit_uppertri_type;
    typedef typename Traits<type>::lowertri_type lowertri_type;
    typedef typename Traits<type>::unit_lowertri_type unit_lowertri_type;
    typedef typename Traits<type>::linearview_type linearview_type;
    typedef typename Traits<type>::realview_type realview_type;
    typedef typename Traits<type>::imagview_type imagview_type;
    typedef typename Traits<type>::nonconj_type nonconj_type;

    typedef typename Traits<type>::reference reference;


    //
    // Constructor
    //

    inline BaseMatrix_Rec_Mutable() {}
    inline BaseMatrix_Rec_Mutable(const BaseMatrix_Rec_Mutable<M>&) {}
    inline ~BaseMatrix_Rec_Mutable() {}


    //
    // Access 
    //

    inline reference operator()(int i, int j)
    {
      CheckRowIndex<mfort>(i,colsize());
      CheckColIndex<mfort>(j,rowsize());
      return ref(i,j);
    }

    // The get_ routines always use CStyle indexing.
    inline row_type get_row(int i) 
    { return row_type(ptr()+i*stepi(),rowsize(),stepj()); }

    inline col_type get_col(int j) 
    { return col_type(ptr()+j*stepj(),colsize(),stepi()); }

    // No need for a get_ routine for diag()
    inline diag_type diag() 
    {
      return diag_type(ptr(),TMV_MIN(colsize(),rowsize()),diagstep());
    }

    inline row_range_type get_row(int i, int j1, int j2) 
    { return row_range_type(ptr()+i*stepi()+j1*stepj(),j2-j1,stepj()); }

    inline col_range_type get_col(int j, int i1, int i2) 
    { return col_range_type(ptr()+j*stepj()+i1*stepi(),i2-i1,stepi()); }

    inline diag_range_type get_diag(int i) 
    {
      return diag_range_type(
          ptr() + (i<0?(-i*stepi()):(i*stepj())),
          (i<0? TMV_MIN(colsize()+i,rowsize()) : TMV_MIN(colsize(),rowsize()-i)),
          diagstep());
    }

    inline diag_range_type get_diag(int i, int j1, int j2) 
    {
      return diag_range_type(
          ptr() + (i<0?(-i*stepi()):(i*stepj())) + j1*diagstep(),
          j2-j1, diagstep());
    }


    // The regular versions respect the indexing style for i and j:
    inline row_type row(int i) 
    {
      CheckRowIndex<mfort>(i,colsize());
      return get_row(i);
    }

    inline row_type operator[](int i) 
    { return row(i); }

    inline col_type col(int j) 
    {
      CheckColIndex<mfort>(j,rowsize());
      return get_col(j);
    }

    inline row_range_type row(int i, int j1, int j2) 
    {
      CheckRowIndex<mfort>(i,colsize());
      CheckColRange<mfort>(j1,j2,rowsize());
      return get_row(i,j1,j2); 
    }

    inline col_range_type col(int j, int i1, int i2) 
    {
      CheckColIndex<mfort>(j,rowsize());
      CheckRowRange<mfort>(i1,i2,colsize());
      return get_col(j,i1,i2); 
    }

    inline diag_range_type diag(int i) 
    {
      CheckDiagIndex<mfort>(i,colsize(),rowsize());
      return get_diag(i);
    }

    inline diag_range_type diag(int i, int j1, int j2) 
    {
      CheckDiagIndex<mfort>(i,j1,j2,colsize(),rowsize());
      return get_diag(i,j1,j2);
    }


    // We need to repeat the const versions so the non-const ones
    // don't clobber them.
    inline value_type operator()(int i, int j) const
    { return base_rec::operator()(i,j); }

    inline const_row_type get_row(int i) const
    { return base_rec::get_row(i); }
    inline const_col_type get_col(int j) const
    { return base_rec::get_col(j); }
    inline const_diag_type diag() const
    { return base_rec::diag(); }
    inline const_row_range_type get_row(int i, int j1, int j2) const
    { return base_rec::get_row(i,j1,j2); }
    inline const_col_range_type get_col(int j, int i1, int i2) const
    { return base_rec::get_col(j,i1,i2); }
    inline const_diag_range_type get_diag(int i) const
    { return base_rec::get_diag(i); }
    inline const_diag_range_type get_diag(int i, int j1, int j2) const
    { return base_rec::get_diag(i,j1,j2); }

    inline const_row_type row(int i) const
    { return base_rec::row(i); }
    inline const_row_type operator[](int i) const
    { return base_rec::row(i); }
    inline const_col_type col(int j) const
    { return base_rec::col(j); }
    inline const_row_range_type row(int i, int j1, int j2) const
    { return base_rec::row(i,j1,j2); }
    inline const_col_range_type col(int j, int i1, int i2) const
    { return base_rec::col(j,i1,i2); }
    inline const_diag_range_type diag(int i) const
    { return base_rec::diag(i); }
    inline const_diag_range_type diag(int i, int j1, int j2) const
    { return base_rec::diag(i,j1,j2); }


    //
    // Op =
    //

    inline type& operator=(BaseMatrix_Rec_Mutable<M>& m2) 
    {
      TMVAssert(colsize() == m2.colsize());
      TMVAssert(rowsize() == m2.rowsize());
      m2.AssignTo(mat());
      return mat(); 
    }

    template <class M2>
    inline type& operator=(const BaseMatrix<M2>& m2) 
    {
      TMVStaticAssert((Sizes<mcolsize,M2::mcolsize>::same));
      TMVStaticAssert((Sizes<mrowsize,M2::mrowsize>::same));
      TMVAssert(colsize() == m2.colsize());
      TMVAssert(rowsize() == m2.rowsize());
      TMVStaticAssert((ShapeTraits2<M2::mshape,mshape>::assignable));
      m2.AssignTo(mat());
      return mat(); 
    }

    inline type& operator=(const value_type x)
    {
      TMVStaticAssert((Sizes<mrowsize,mcolsize>::same));
      TMVAssert(colsize() == rowsize());
      SetToIdentity(x);
      return mat();
    }

    inline ListAssigner<value_type,typename linearview_type::iterator> 
    operator<<(value_type x)
    {
      TMVAssert(this->CanLinearize());
      linearview_type v = LinearView();
      typedef typename linearview_type::iterator iterator;
      return (ListAssigner<value_type,iterator>(v.begin(),v.size()),x); 
    }


    //
    // Modifying Functions
    //

    // First the ones from BaseMatrix_Mutable:
    inline type& Zero()
    { tmv::Zero(mat()); return mat(); }

    inline type& SetAllTo(value_type val)
    { tmv::SetAllTo(mat(),val); return mat(); }

    inline type& AddToAll(value_type val)
    { tmv::AddToAll(mat(),val); return mat(); }

    inline type& Clip(real_type thresh)
    { tmv::Clip(mat(),thresh); return mat(); }

    template <class F>
    inline type& ApplyToAll(const F& f)
    { tmv::ApplyToAll(mat(),f); return mat(); }

    inline type& ConjugateSelf()
    { tmv::ConjugateSelf(mat()); return mat(); }

    // Some more that are added for Rec shape:
    inline type& TransposeSelf() 
    { tmv::TransposeSelf(mat()); return mat(); }

    inline type& SetToIdentity(const value_type x=value_type(1))
    {
      TMVStaticAssert((Sizes<mrowsize,mcolsize>::same));
      TMVAssert(colsize() == rowsize());
      this->Zero(); diag().SetAllTo(x);
      return mat();
    }

    inline type& CSwapRows(int i1, int i2) 
    {
      if (i1 != i2) {
        tmv::Swap(get_row(i1),get_row(i2));
      }
      return mat();
    }
    inline type& SwapRows(int i1, int i2) 
    {
      CheckRowIndex<mfort>(i1,colsize());
      CheckRowIndex<mfort>(i2,colsize());
      return CSwapRows(i1,i2);
    }

    inline type& CSwapCols(int j1, int j2) 
    {
      if (j1 != j2) {
        tmv::Swap(get_col(j1),get_col(j2));
      }
      return mat();
    }
    inline type& SwapCols(int j1, int j2) 
    {
      CheckColIndex<mfort>(j1,rowsize());
      CheckColIndex<mfort>(j2,rowsize());
      return CSwapCols(j1,j2);
    }

    inline type& CPermuteRows(const int*const p, int i1, int i2) 
    { tmv::PermuteRows(mat(),p,i1,i2); return mat(); }
    inline type& PermuteRows(const int*const p, int i1, int i2) 
    {
      CheckRowRange<mfort>(i1,i2,colsize());
      return CPermuteRows(p,i1,i2);
    }
    inline type& PermuteRows(const int*const p) 
    { CPermuteRows(p,0,colsize()); return mat(); }

    inline type& CReversePermuteRows(const int*const p, int i1, int i2) 
    { tmv::ReversePermuteRows(mat(),p,i1,i2); }
    inline type& ReversePermuteRows(const int*const p, int i1, int i2) 
    {
      CheckRowRange<mfort>(i1,i2,colsize());
      return CReversePermuteRows(p,i1,i2);
    }
    inline type& ReversePermuteRows(const int*const p) 
    { CReversePermuteRows(p,0,colsize()); return mat(); }

    inline type& PermuteCols(const int*const p, int j1, int j2) 
    {
      CheckColRange<mfort>(j1,j2,rowsize());
      Transpose().CPermuteRows(p,j1,j2);
      return mat();
    }
    inline type& PermuteCols(const int*const p) 
    { Transpose().CPermuteRows(p,0,rowsize()); return mat(); }

    inline type& ReversePermuteCols(const int*const p, int j1, int j2) 
    {
      CheckColRange<mfort>(j1,j2,rowsize());
      Transpose().CReversePermuteRows(p,j1,j2);
      return mat();
    }
    inline type& ReversePermuteCols(const int*const p) 
    { Transpose().CReversePermuteRows(p,0,rowsize()); return mat(); }


    //
    // SubMatrix, etc.
    //

    // These versions always uses CStyle
    inline submatrix_type CSubMatrix(int i1, int i2, int j1, int j2) 
    {
      return submatrix_type(
          ptr()+i1*stepi()+j1*stepj(), i2-i1, j2-j1, stepi(), stepj()); 
    }

    inline submatrix_step_type CSubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) 
    {
      return submatrix_step_type(
          ptr()+i1*stepi()+j1*stepj(), (i2-i1)/istep, (j2-j1)/jstep,
          istep*stepi(), jstep*stepj());
    }

    inline subvector_type CSubVector(
        int i, int j, int istep, int jstep, int s) 
    {
      return subvector_type(
          ptr()+i*stepi()+j*stepj(), s, istep*stepi() + jstep*stepj());
    }

    inline colpair_type CColPair(int j1, int j2) 
    {
      return colpair_type(
          ptr()+j1*stepj(), colsize(), 2, stepi(), (j2-j1)*stepj());
    }

    inline rowpair_type CRowPair(int i1, int i2) 
    {
      return rowpair_type(
          ptr()+i1*stepi(), 2, rowsize(), (i2-i1)*stepi(), stepj());
    }

    inline cols_type CCols(int j1, int j2) 
    {
      return cols_type(
          ptr()+j1*stepj(), colsize(), j2-j1, stepi(), stepj());
    }

    inline rows_type CRows(int i1, int i2) 
    {
      return rows_type(
          ptr()+i1*stepi(), i2-i1, rowsize(), stepi(), stepj());
    }


    // These check the indices according the the indexing style being
    // used, and then calls the above CStyle versions.
    inline submatrix_type SubMatrix(
        int i1, int i2, int j1, int j2) 
    {
      CheckRowRange<mfort>(i1,i2,colsize());
      CheckColRange<mfort>(j1,j2,rowsize());
      return CSubMatrix(i1,i2,j1,j2);
    }

    inline submatrix_step_type SubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) 
    {
      CheckRowRange<mfort>(i1,i2,istep,colsize());
      CheckColRange<mfort>(j1,j2,jstep,rowsize());
      return CSubMatrix(i1,i2,j1,j2,istep,jstep);
    }

    inline subvector_type SubVector(
        int i, int j, int istep, int jstep, int s) 
    {
      CheckMatSubVector<mfort>(i,j,istep,jstep,s,colsize(),rowsize());
      return CSubVector(i,j,istep,jstep,s); 
    }

    inline colpair_type ColPair(int j1, int j2) 
    {
      CheckColIndex<mfort>(j1,rowsize());
      CheckColIndex<mfort>(j2,rowsize());
      return CColPair(j1,j2);
    }

    inline rowpair_type RowPair(int i1, int i2) 
    {
      CheckRowIndex<mfort>(i1,colsize());
      CheckRowIndex<mfort>(i2,colsize());
      return CRowPair(i1,i2);
    }

    inline cols_type Cols(int j1, int j2) 
    {
      CheckColRange<mfort>(j1,j2,rowsize());
      return CCols(j1,j2);
    }

    inline rows_type Rows(int i1, int i2) 
    {
      CheckRowRange<mfort>(i1,i2,colsize());
      return CRows(i1,i2);
    }


    // Repeat the const versions:
    inline const_submatrix_type CSubMatrix(
        int i1, int i2, int j1, int j2) const
    { return base_rec::CSubMatrix(i1,i2,j1,j2); }
    inline const_submatrix_step_type CSubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) const
    { return base_rec::CSubMatrix(i1,i2,j1,j2,istep,jstep); }
    inline const_subvector_type CSubVector(
        int i, int j, int istep, int jstep, int s) const
    { return base_rec::CSubVector(i,j,istep,jstep,s); }
    inline const_colpair_type CColPair(int j1, int j2) const
    { return base_rec::CColPair(j1,j2); }
    inline const_rowpair_type CRowPair(int i1, int i2) const
    { return base_rec::CRowPair(i1,i2); }
    inline const_cols_type CCols(int j1, int j2) const
    { return base_rec::CCols(j1,j2); }
    inline const_rows_type CRows(int i1, int i2) const
    { return base_rec::CRows(i1,i2); }

    inline const_submatrix_type SubMatrix(
        int i1, int i2, int j1, int j2) const
    { return base_rec::SubMatrix(i1,i2,j1,j2); }
    inline const_submatrix_step_type SubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) const
    { return base_rec::SubMatrix(i1,i2,j1,j2,istep,jstep); }
    inline const_subvector_type SubVector(
        int i, int j, int istep, int jstep, int s) const
    { return base_rec::SubVector(i,j,istep,jstep,s); }
    inline const_colpair_type ColPair(int j1, int j2) const
    { return base_rec::ColPair(j1,j2); }
    inline const_rowpair_type RowPair(int i1, int i2) const
    { return base_rec::RowPair(i1,i2); }
    inline const_cols_type Cols(int j1, int j2) const
    { return base_rec::Cols(j1,j2); }
    inline const_rows_type Rows(int i1, int i2) const
    { return base_rec::Rows(i1,i2); }


    //
    // Views
    //

    inline view_type View() 
    { return view_type(ptr(),colsize(),rowsize(),stepi(),stepj()); }

    inline cview_type CView() 
    { return View(); }

    inline fview_type FView() 
    { return View(); }

    inline xview_type XView() 
    { return View(); }

    inline cmview_type CMView() 
    { return View(); }

    inline rmview_type RMView() 
    { return View(); }

    inline transpose_type Transpose() 
    { return transpose_type(ptr(),rowsize(),colsize(),stepj(),stepi()); }

    inline conjugate_type Conjugate() 
    { return conjugate_type(ptr(),colsize(),rowsize(),stepi(),stepj()); }

    inline adjoint_type Adjoint() 
    { return adjoint_type(ptr(),rowsize(),colsize(),stepj(),stepi()); }

    inline uppertri_type UpperTri() 
    { return uppertri_type(ptr(),rowsize(),stepi(),stepj()); }

    inline unit_uppertri_type UnitUpperTri() 
    { return unit_uppertri_type(ptr(),rowsize(),stepi(),stepj()); }

    inline lowertri_type LowerTri() 
    { return lowertri_type(ptr(),colsize(),stepi(),stepj()); }

    inline unit_lowertri_type UnitLowerTri() 
    { return unit_lowertri_type(ptr(),colsize(),stepi(),stepj()); }

    inline linearview_type LinearView() 
    {
      TMVAssert(this->CanLinearize());
      return linearview_type(ptr(),ls(),1);
    }

    inline realview_type Real() 
    {
      return realview_type(
          reinterpret_cast<real_type*>(ptr()), colsize(), rowsize(),
          misreal ? stepi() : 2*stepi(), misreal ? stepj() : 2*stepj());
    }

    inline imagview_type Imag() 
    {
      TMVStaticAssert(miscomplex);
      return imagview_type(
          reinterpret_cast<real_type*>(ptr())+1, colsize(), rowsize(),
          misreal ? stepi() : 2*stepi(), misreal ? stepj() : 2*stepj());
    }

    inline nonconj_type NonConj()
    { return nonconj_type(ptr(),colsize(),rowsize(),stepi(),stepj()); }


    // Repeat the const versions:
    inline const_view_type View() const
    { return base_rec::View(); }
    inline const_cview_type CView() const
    { return base_rec::View(); }
    inline const_fview_type FView() const
    { return base_rec::View(); }
    inline const_xview_type XView() const
    { return base_rec::View(); }
    inline const_cmview_type CMView() const
    { return base_rec::View(); }
    inline const_rmview_type RMView() const
    { return base_rec::View(); }
    inline const_transpose_type Transpose() const
    { return base_rec::Transpose(); }
    inline const_conjugate_type Conjugate() const
    { return base_rec::Conjugate(); }
    inline const_adjoint_type Adjoint() const
    { return base_rec::Adjoint(); }
    inline const_uppertri_type UpperTri() const
    { return base_rec::UpperTri(); }
    inline const_unit_uppertri_type UnitUpperTri() const
    { return base_rec::UnitUpperTri(); }
    inline const_lowertri_type LowerTri() const
    { return base_rec::LowerTri(); }
    inline const_unit_lowertri_type UnitLowerTri() const
    { return base_rec::UnitLowerTri(); }
    inline const_linearview_type LinearView() const
    { return base_rec::LinearView(); }
    inline const_realview_type Real() const
    { return base_rec::Real(); }
    inline const_imagview_type Imag() const
    { return base_rec::Imag(); }
    inline const_nonconj_type NonConj() const
    { return base_rec::NonConj(); }



    //
    // I/O
    //

    inline void Read(std::istream& is)
    {
      cview_type mcv = CView();
      tmv::Read(is,mcv); 
    }

    //
    // Auxilliary routines
    //

    inline const type& mat() const
    { return *static_cast<const type*>(this); }
    inline type& mat()
    { return *static_cast<type*>(this); }

    inline bool isconj() const { return mconj; }
    inline int diagstep() const 
    { return mdiagstep == UNKNOWN ? stepi() + stepj() : mdiagstep; }
    inline bool CanLinearize() const 
    { return (AuxCanLinearize<mcanlin,mcolmajor,mrowmajor,M>::ok(mat())); }

    // Note that these last functions need to be defined in a more derived
    // class than this, or an infinite loop will result when compiling.
    // Also, cref and cptr from above.

    inline size_t colsize() const { return mat().colsize(); }
    inline size_t rowsize() const { return mat().rowsize(); }
    inline size_t ls() const { return mat().ls(); }
    inline int stepi() const { return mat().stepi(); }
    inline int stepj() const { return mat().stepj(); }
    inline bool isrm() const { return mat().isrm(); }
    inline bool iscm() const { return mat().iscm(); }

    inline value_type* ptr() { return mat().ptr(); }
    inline reference ref(int i, int j) { return mat().ref(i,j); }

  }; // BaseMatrix_Rec_Mutable

  //
  // Zero
  //

  template <int algo, class M> struct ZeroM_Helper;

  // algo 1: Linearize to vector version
  template <class M> 
  struct ZeroM_Helper<1,M> 
  {
    static inline void call(M& m) 
    { m.LinearView().Zero(); } 
  };

  // algo 2: RowMajor
  template <class M> 
  struct ZeroM_Helper<2,M> 
  {
    static inline void call(M& m1) 
    {
      const int m = m1.colsize();
      for(int i=0;i<m;++i) m1.get_row(i).Zero();
    } 
  };

  // algo 3: ColMajor
  template <class M> 
  struct ZeroM_Helper<3,M> 
  {
    static inline void call(M& m) 
    {
      const int n = m.rowsize();
      for(int j=0;j<n;++j) m.get_col(j).Zero();
    } 
  };

  // algo 4: Unknown sizes, determine which algorithm to use
  template <class M> 
  struct ZeroM_Helper<4,M> 
  {
    static inline void call(M& m) 
    {
      enum { algo2 = M::mcolmajor ? 3 : 2 };
#if TMV_OPT >= 2
      if (m.CanLinearize())
        ZeroM_Helper<1,M>::call(m);
      else
#endif
        ZeroM_Helper<algo2,M>::call(m);
    } 
  };

  template <class M>
  inline void Zero(BaseMatrix_Rec_Mutable<M>& m)
  {
#if TMV_OPT == 0
    enum { algo = M::mcolmajor ? 3 : 2 };
#else
    enum { algo = M::mcanlin ? 1 : 4 };
#endif
    ZeroM_Helper<algo,M>::call(m.mat());
  }

  //
  // SetAllTo
  //

  template <int algo, class M, class T> struct SetAllToM_Helper;

  // algo 1: Linearize to vector version
  template <class M, class T> 
  struct SetAllToM_Helper<1,M,T> // algo 1, linearize
  { 
    static inline void call(M& m, const T& val) 
    { m.LinearView().SetAllTo(val); } 
  };

  // algo 2: RowMajor
  template <class M, class T> 
  struct SetAllToM_Helper<2,M,T>
  {
    static void call(M& m1, const T& val) 
    {
      const int m = m1.colsize();
      for(int i=0;i<m;++i) m1.get_row(i).SetAllTo(val);
    } 
  };

  // algo 3: ColMajor
  template <class M, class T> 
  struct SetAllToM_Helper<3,M,T> 
  {
    static void call(M& m, const T& val) 
    {
      const int n = m.rowsize();
      for(int j=0;j<n;++j) m.get_col(j).SetAllTo(val);
    } 
  };

  // algo 4: Unknown sizes, determine which algorithm to use
  template <class M, class T> 
  struct SetAllToM_Helper<4,M,T> 
  {
    static void call(M& m, const T& val) 
    {
      enum { algo2 = M::mcolmajor ? 3 : 2 };
#if TMV_OPT >= 2
      if (m.CanLinearize())
        SetAllToM_Helper<1,M,T>::call(m,val);
      else
#endif
        SetAllToM_Helper<algo2,M,T>::call(m,val);
    } 
  };

  template <class M, class T>
  inline void SetAllTo(BaseMatrix_Rec_Mutable<M>& m, const T& val)
  {
#if TMV_OPT == 0
    enum { algo = M::mcolmajor ? 3 : 2 };
#else
    enum { algo = M::mcanlin ? 1 : 4 };
#endif
    SetAllToM_Helper<algo,M,T>::call(m.mat(),val);
  }

  //
  // AddToAll
  //

  template <int algo, class M, class T> struct AddToAllM_Helper;

  // algo 1: Linearize to vector version
  template <class M, class T> 
  struct AddToAllM_Helper<1,M,T> // algo 1, linearize
  { 
    static inline void call(M& m, const T& val) 
    { m.LinearView().AddToAll(val); } 
  };

  // algo 2: RowMajor
  template <class M, class T> 
  struct AddToAllM_Helper<2,M,T> 
  {
    static void call(M& m1, const T& val) 
    {
      const int m = m1.colsize();
      for(int i=0;i<m;++i) m1.get_row(i).AddToAll(val);
    } 
  };

  // algo 3: ColMajor
  template <class M, class T> 
  struct AddToAllM_Helper<3,M,T> 
  {
    static void call(M& m, const T& val) 
    {
      const int n = m.rowsize();
      for(int j=0;j<n;++j) m.get_col(j).AddToAll(val);
    } 
  };

  // algo 4: Unknown sizes, determine which algorithm to use
  template <class M, class T> 
  struct AddToAllM_Helper<4,M,T> 
  {
    static void call(M& m, const T& val) 
    {
      enum { algo2 = M::mcolmajor ? 3 : 2 };
#if TMV_OPT >= 2
      if (m.CanLinearize())
        AddToAllM_Helper<1,M,T>::call(m,val);
      else
#endif
        AddToAllM_Helper<algo2,M,T>::call(m,val);
    } 
  };

  template <class M, class T>
  inline void AddToAll(BaseMatrix_Rec_Mutable<M>& m, const T& val)
  {
#if TMV_OPT == 0
    enum { algo = M::mcolmajor ? 3 : 2 };
#else
    enum { algo = M::mcanlin ? 1 : 4 };
#endif
    AddToAllM_Helper<algo,M,T>::call(m.mat(),val);
  }

  //
  // Clip
  //

  template <int algo, class M, class RT> struct ClipM_Helper;

  // algo 1: Linearize to vector version
  template <class M, class RT> 
  struct ClipM_Helper<1,M,RT> // algo 1, linearize
  {
    static inline void call(M& m, const RT& thresh) 
    { m.LinearView().Clip(thresh); } 
  };

  // algo 2: RowMajor
  template <class M, class RT> 
  struct ClipM_Helper<2,M,RT> 
  {
    static void call(M& m1, const RT& thresh) 
    {
      const int m = m1.colsize();
      for(int i=0;i<m;++i) m1.get_row(i).Clip(thresh);
    } 
  };

  // algo 3: ColMajor
  template <class M, class RT> 
  struct ClipM_Helper<3,M,RT> 
  {
    static void call(M& m, const RT& thresh) 
    {
      const int n = m.rowsize();
      for(int j=0;j<n;++j) m.get_col(j).Clip(thresh);
    } 
  };

  // algo 4: Unknown sizes, determine which algorithm to use
  template <class M, class RT> 
  struct ClipM_Helper<4,M,RT> 
  {
    static void call(M& m, const RT& thresh) 
    {
      enum { algo2 = M::mcolmajor ? 3 : 2 };
#if TMV_OPT >= 2
      if (m.CanLinearize())
        ClipM_Helper<1,M,RT>::call(m,thresh);
      else
#endif
        ClipM_Helper<algo2,M,RT>::call(m,thresh);
    } 
  };

  template <class M, class RT>
  inline void Clip(BaseMatrix_Rec_Mutable<M>& m, const RT& thresh)
  {
#if TMV_OPT == 0
    enum { algo = M::mcolmajor ? 3 : 2 };
#else
    enum { algo = M::mcanlin ? 1 : 4 };
#endif
    ClipM_Helper<algo,M,RT>::call(m.mat(),thresh);
  }

  //
  // ApplyToAll
  //

  template <int algo, class M, class F> struct ApplyToAllM_Helper;

  // algo 1: Linearize to vector version
  template <class M, class F> 
  struct ApplyToAllM_Helper<1,M,F> // algo 1, linearize
  { 
    static inline void call(M& m, const F& f) 
    { m.LinearView().ApplyToAll(f); } 
  };

  // algo 2: RowMajor
  template <class M, class F> 
  struct ApplyToAllM_Helper<2,M,F> 
  {
    static void call(M& m1, const F& f) 
    {
      const int m = m1.colsize();
      for(int i=0;i<m;++i) m1.get_row(i).ApplyToAll(f);
    } 
  };

  // algo 3: ColMajor
  template <class M, class F> 
  struct ApplyToAllM_Helper<3,M,F> 
  {
    static void call(M& m, const F& f) 
    {
      const int n = m.rowsize();
      for(int j=0;j<n;++j) m.get_col(j).ApplyToAll(f);
    } 
  };

  // algo 4: Unknown sizes, determine which algorithm to use
  template <class M, class F> 
  struct ApplyToAllM_Helper<4,M,F>
  {
    static void call(M& m, const F& f) 
    {
      enum { algo2 = M::mcolmajor ? 3 : 2 };
#if TMV_OPT >= 2
      if (m.CanLinearize())
        ApplyToAllM_Helper<1,M,F>::call(m,f);
      else
#endif
        ApplyToAllM_Helper<algo2,M,F>::call(m,f);
    } 
  };

  template <class M, class F>
  inline void ApplyToAll(BaseMatrix_Mutable<M>& m, const F& f)
  {
#if TMV_OPT == 0
    enum { algo = M::mcolmajor ? 3 : 2 };
#else
    enum { algo = M::mcanlin ? 1 : 4 };
#endif
    ApplyToAllM_Helper<algo,M,F>::call(m.mat(),f);
  }

  //
  // ConjugateSelf
  //

  template <int algo, class M> struct ConjugateM_Helper;

  // algo 0: Not complex, nothing to do
  template <class M>
  struct ConjugateM_Helper<0,M>
  { static inline void call(M& ) {} };

  // algo 1: Linearize to vector version
  template <class M>
  struct ConjugateM_Helper<1,M>
  {
    static inline void call(M& m)
    { m.LinearView().ConjugateSelf(); }
  };

  // In TMV_ProdXM.h
  template <int ix, class T, class M>
  inline void MultXM(const Scaling<ix,T>& x, BaseMatrix_Mutable<M>& m);

  // algo 2: m.Imag() *= -1
  template <class M>
  struct ConjugateM_Helper<2,M>
  {
    static inline void call(M& m)
    {
      typedef typename M::real_type RT;
      typedef typename M::imagview_type Mi;
      const Scaling<-1,RT> mone;
      Mi mi = m.Imag();
      MultXM(mone,mi);
    }
  };

  // algo 4: Unknown sizes, determine which algorithm to use
  template <class M>
  struct ConjugateM_Helper<4,M>
  {
    static void call(M& m)
    {
#if TMV_OPT >= 2
      if (m.CanLinearize())
        ConjugateM_Helper<1,M>::call(m);
      else
#endif
        ConjugateM_Helper<2,M>::call(m);
    }
  };

  template <class M>
  inline void ConjugateSelf(BaseMatrix_Rec_Mutable<M>& m)
  {
#if TMV_OPT == 0
    enum { algo = M::misreal ? 0 : 2 };
#else
    enum { algo = M::misreal ? 0 : M::mcanlin ? 1 : 4 };
#endif
    ConjugateM_Helper<algo,M>::call(m.mat());
  }


  //
  // TypeText 
  //

  template <class M>
  static std::string TypeText(const BaseMatrix_Rec<M>& m)
  {
    std::ostringstream s;
    s << "BaseMatrix_Rec< "<<TypeText(m.mat())<<" >";
    return s.str();
  }

  template <class M>
  static std::string TypeText(const BaseMatrix_Rec_Mutable<M>& m)
  {
    std::ostringstream s;
    s << "BaseMatrix_Rec_Mutable< "<<TypeText(m.mat())<<" >";
    return s.str();
  }

} // namespace tmv

#endif
