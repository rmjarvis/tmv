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
// This file defines the BaseMatrix_Tri and BaseMatrix_Tri classes.
//
// See TMV_TriMatrix.h for the functions that are defined for these objects.
//

#ifndef TMV_BaseMatrix_Tri_H
#define TMV_BaseMatrix_Tri_H

#include "tmv/TMV_BaseMatrix.h"
#include "tmv/TMV_MultVV.h"
#include "tmv/TMV_BaseMatrix_Rec.h"
#include "tmv/TMV_BaseMatrix_Diag.h"

namespace tmv {

  template <class M>
  class BaseMatrix_Tri;
  template <class M>
  class BaseMatrix_Tri_Mutable;

  template <class T, DiagType D=NonUnitDiag, StorageType S=ColMajor, IndexStyle I=CStyle>
  class UpperTriMatrix;
  template <class T, DiagType D=NonUnitDiag, int Si=UNKNOWN, int Sj=UNKNOWN, bool C=false, IndexStyle I=CStyle>
  class ConstUpperTriMatrixView;
  template <class T, DiagType D=NonUnitDiag, int Si=UNKNOWN, int Sj=UNKNOWN, bool C=false, IndexStyle I=CStyle>
  class UpperTriMatrixView;
  template <class T, int N, DiagType D=NonUnitDiag, StorageType S=ColMajor, IndexStyle I=CStyle>
  class SmallUpperTriMatrix;
  template <class T, int N, DiagType D, int Si, int Sj, bool C=false, IndexStyle I=CStyle>
  class ConstSmallUpperTriMatrixView;
  template <class T, int N, DiagType D, int Si, int Sj, bool C=false, IndexStyle I=CStyle>
  class SmallUpperTriMatrixView;

  template <class T, DiagType D=NonUnitDiag, StorageType S=ColMajor, IndexStyle I=CStyle>
  class LowerTriMatrix;
  template <class T, DiagType D=NonUnitDiag, int Si=UNKNOWN, int Sj=UNKNOWN, bool C=false, IndexStyle I=CStyle>
  class ConstLowerTriMatrixView;
  template <class T, DiagType D=NonUnitDiag, int Si=UNKNOWN, int Sj=UNKNOWN, bool C=false, IndexStyle I=CStyle>
  class LowerTriMatrixView;
  template <class T, int N, DiagType D=NonUnitDiag, StorageType S=ColMajor, IndexStyle I=CStyle>
  class SmallLowerTriMatrix;
  template <class T, int N, DiagType D, int Si, int Sj, bool C=false, IndexStyle I=CStyle>
  class ConstSmallLowerTriMatrixView;
  template <class T, int N, DiagType D, int Si, int Sj, bool C=false, IndexStyle I=CStyle>
  class SmallLowerTriMatrixView;

  template <class T, DiagType D=NonUnitDiag, StorageType S=ColMajor>
  class UpperTriMatrixF;
  template <class T, DiagType D=NonUnitDiag, int Si=UNKNOWN, int Sj=UNKNOWN, bool C=false>
  class ConstUpperTriMatrixViewF;
  template <class T, DiagType D=NonUnitDiag, int Si=UNKNOWN, int Sj=UNKNOWN, bool C=false>
  class UpperTriMatrixViewF;
  template <class T, int N, DiagType D=NonUnitDiag, StorageType S=ColMajor>
  class SmallUpperTriMatrixF;
  template <class T, int N, DiagType D, int Si, int Sj, bool C=false>
  class ConstSmallUpperTriMatrixViewF;
  template <class T, int N, DiagType D, int Si, int Sj, bool C=false>
  class SmallUpperTriMatrixViewF;

  template <class T, DiagType D=NonUnitDiag, StorageType S=ColMajor>
  class LowerTriMatrixF;
  template <class T, DiagType D=NonUnitDiag, int Si=UNKNOWN, int Sj=UNKNOWN, bool C=false>
  class ConstLowerTriMatrixViewF;
  template <class T, DiagType D=NonUnitDiag, int Si=UNKNOWN, int Sj=UNKNOWN, bool C=false>
  class LowerTriMatrixViewF;
  template <class T, int N, DiagType D=NonUnitDiag, StorageType S=ColMajor>
  class SmallLowerTriMatrixF;
  template <class T, int N, DiagType D, int Si, int Sj, bool C=false>
  class ConstSmallLowerTriMatrixViewF;
  template <class T, int N, DiagType D, int Si, int Sj, bool C=false>
  class SmallLowerTriMatrixViewF;


  // Specify ExactSameStorage for triangle matrices:
  template <class M>
  inline bool ExactSameStorage(
      const BaseMatrix_Tri<M>& m1, const BaseMatrix_Tri<M>& m2)
  { return m1.stepi() == m2.stepi() && m1.stepj() == m2.stepj(); }

  template <int shape>
  inline void CheckTri(int i, int j)
  { // UpperTri
    TMVAssert(i <= j && "i,j must be in upper triangle");
  }
  template <>
  inline void CheckTri<LowerTri>(int i, int j)
  { // LowerTri
    TMVAssert(i >= j && "i,j must be in lower triangle");
  }
  template <>
  inline void CheckTri<UnitUpperTri>(int i, int j)
  { // UnitUpperTri
    TMVAssert(i < j && "i,j must be in strict upper triangle");
  }
  template <>
  inline void CheckTri<UnitLowerTri>(int i, int j)
  { // UnitLowerTri
    TMVAssert(i > j && "i,j must be in strict upper triangle");
  }
  template <int shape>
  inline void CheckDiagTri(int i)
  { // UpperTri
    TMVAssert(i >= 0 && 
        "sub-diagonal not allowed for upper triangle matrix");
  }
  template <>
  inline void CheckDiagTri<LowerTri>(int i)
  { // LowerTri
    TMVAssert(i <= 0 && 
        "super-diagonal not allowed for lower triangle matrix");
  }
  template <>
  inline void CheckDiagTri<UnitUpperTri>(int i)
  { // UpperTri
    TMVAssert(i >= 0 && 
        "sub-diagonal not allowed for upper triangle matrix");
    TMVAssert(i != 0 && 
        "main diagonal not allowed for UnitDiag triangle matrix");
  }
  template <>
  inline void CheckDiagTri<UnitLowerTri>(int i)
  { // LowerTri
    TMVAssert(i <= 0 && 
        "super-diagonal not allowed for lower triangle matrix");
    TMVAssert(i != 0 && 
        "main diagonal not allowed for UnitDiag triangle matrix");
  }

  template <class T, int cs, int rs, bool rm, bool fort>
  struct MCopyHelper<T,UpperTri,cs,rs,rm,fort>
  { typedef SmallUpperTriMatrix<T,cs,NonUnitDiag,rm?RowMajor:ColMajor,fort?FortranStyle:CStyle> type; };
  template <class T, int rs, bool rm, bool fort>
  struct MCopyHelper<T,UpperTri,UNKNOWN,rs,rm,fort>
  { typedef SmallUpperTriMatrix<T,rs,NonUnitDiag,rm?RowMajor:ColMajor,fort?FortranStyle:CStyle> type; };
  template <class T, int cs, bool rm, bool fort>
  struct MCopyHelper<T,UpperTri,cs,UNKNOWN,rm,fort>
  { typedef SmallUpperTriMatrix<T,cs,NonUnitDiag,rm?RowMajor:ColMajor,fort?FortranStyle:CStyle> type; };
  template <class T, bool rm, bool fort>
  struct MCopyHelper<T,UpperTri,UNKNOWN,UNKNOWN,rm,fort>
  { typedef UpperTriMatrix<T,NonUnitDiag,rm?RowMajor:ColMajor,fort?FortranStyle:CStyle> type; };

  template <class T, int cs, int rs, bool rm, bool fort>
  struct MCopyHelper<T,LowerTri,cs,rs,rm,fort>
  { typedef SmallLowerTriMatrix<T,cs,NonUnitDiag,rm?RowMajor:ColMajor,fort?FortranStyle:CStyle> type; };
  template <class T, int rs, bool rm, bool fort>
  struct MCopyHelper<T,LowerTri,UNKNOWN,rs,rm,fort>
  { typedef SmallLowerTriMatrix<T,rs,NonUnitDiag,rm?RowMajor:ColMajor,fort?FortranStyle:CStyle> type; };
  template <class T, int cs, bool rm, bool fort>
  struct MCopyHelper<T,LowerTri,cs,UNKNOWN,rm,fort>
  { typedef SmallLowerTriMatrix<T,cs,NonUnitDiag,rm?RowMajor:ColMajor,fort?FortranStyle:CStyle> type; };
  template <class T, bool rm, bool fort>
  struct MCopyHelper<T,LowerTri,UNKNOWN,UNKNOWN,rm,fort>
  { typedef LowerTriMatrix<T,NonUnitDiag,rm?RowMajor:ColMajor,fort?FortranStyle:CStyle> type; };

  template <class T, int cs, int rs, bool rm, bool fort>
  struct MCopyHelper<T,UnitUpperTri,cs,rs,rm,fort>
  { typedef SmallUpperTriMatrix<T,cs,UnitDiag,rm?RowMajor:ColMajor,fort?FortranStyle:CStyle> type; };
  template <class T, int rs, bool rm, bool fort>
  struct MCopyHelper<T,UnitUpperTri,UNKNOWN,rs,rm,fort>
  { typedef SmallUpperTriMatrix<T,rs,UnitDiag,rm?RowMajor:ColMajor,fort?FortranStyle:CStyle> type; };
  template <class T, int cs, bool rm, bool fort>
  struct MCopyHelper<T,UnitUpperTri,cs,UNKNOWN,rm,fort>
  { typedef SmallUpperTriMatrix<T,cs,UnitDiag,rm?RowMajor:ColMajor,fort?FortranStyle:CStyle> type; };
  template <class T, bool rm, bool fort>
  struct MCopyHelper<T,UnitUpperTri,UNKNOWN,UNKNOWN,rm,fort>
  { typedef UpperTriMatrix<T,UnitDiag,rm?RowMajor:ColMajor,fort?FortranStyle:CStyle> type; };

  template <class T, int cs, int rs, bool rm, bool fort>
  struct MCopyHelper<T,UnitLowerTri,cs,rs,rm,fort>
  { typedef SmallLowerTriMatrix<T,cs,UnitDiag,rm?RowMajor:ColMajor,fort?FortranStyle:CStyle> type; };
  template <class T, int rs, bool rm, bool fort>
  struct MCopyHelper<T,UnitLowerTri,UNKNOWN,rs,rm,fort>
  { typedef SmallLowerTriMatrix<T,rs,UnitDiag,rm?RowMajor:ColMajor,fort?FortranStyle:CStyle> type; };
  template <class T, int cs, bool rm, bool fort>
  struct MCopyHelper<T,UnitLowerTri,cs,UNKNOWN,rm,fort>
  { typedef SmallLowerTriMatrix<T,cs,UnitDiag,rm?RowMajor:ColMajor,fort?FortranStyle:CStyle> type; };
  template <class T, bool rm, bool fort>
  struct MCopyHelper<T,UnitLowerTri,UNKNOWN,UNKNOWN,rm,fort>
  { typedef LowerTriMatrix<T,UnitDiag,rm?RowMajor:ColMajor,fort?FortranStyle:CStyle> type; };


  // Defined in TMV_CopyU.h
  template <class M1, class M2> 
  inline void Copy(
      const BaseMatrix_Tri<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2);
  template <class M1, class M2> 
  inline void Copy(
      const BaseMatrix_Tri<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2);

  // Defined in TMV_SwapU.h
  template <class M1, class M2>
  inline void Swap(
      BaseMatrix_Tri_Mutable<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2);

  // Defined in TMV_NormU.h
  template <class M>
  inline typename M::real_type NormSq(const BaseMatrix_Tri<M>& m);
  template <class M>
  inline typename M::real_type NormSq(const BaseMatrix_Tri<M>& m,
      const typename M::real_type scale);
  template <class M>
  inline typename M::real_type NormF(const BaseMatrix_Tri<M>& m);
  template <class M>
  inline typename M::real_type Norm1(const BaseMatrix_Tri<M>& m);
  template <class M>
  inline typename M::real_type NormInf(const BaseMatrix_Tri<M>& m);
  template <class M>
  inline typename M::value_type SumElements(const BaseMatrix_Tri<M>& m);
  template <class M>
  inline typename M::real_type SumAbsElements(const BaseMatrix_Tri<M>& m);
  template <class M>
  inline typename M::real_type MaxAbsElement(const BaseMatrix_Tri<M>& m);

  // Defined in TMV_TriMatrixIO.h
  template <class M>
  inline void WriteCompact(std::ostream& os, const BaseMatrix_Tri<M>& m);
  template <class M>
  inline void WriteCompact(std::ostream& os, const BaseMatrix_Tri<M>& m,
      typename M::real_type thresh) ;
  template <class M>
  inline void Read(std::istream& is, BaseMatrix_Tri_Mutable<M>& m);

  // Defined below:
  template <class M>
  inline void Zero(BaseMatrix_Tri_Mutable<M>& m);
  template <class M, class T>
  inline void SetAllTo(BaseMatrix_Tri_Mutable<M>& m, const T& val);
  template <class M, class T>
  inline void AddToAll(BaseMatrix_Tri_Mutable<M>& m, const T& val);
  template <class M, class RT>
  inline void Clip(BaseMatrix_Tri_Mutable<M>& m, const RT& thresh);
  template <class M, class F>
  inline void ApplyToAll(BaseMatrix_Tri_Mutable<M>& m, const F& f);
  template <class M>
  inline void ConjugateSelf(BaseMatrix_Tri_Mutable<M>& m);
  template <class M>
  inline void InvertTri(BaseMatrix_Tri_Mutable<M>& m);


  template <class M>
  class BaseMatrix_Tri :
    public BaseMatrix_Calc<M>
  {
  public:
    typedef M type;

    typedef typename Traits<type>::value_type value_type;
    typedef typename Traits<type>::calc_type calc_type;
    typedef typename Traits<type>::eval_type eval_type;
    typedef typename Traits<type>::copy_type copy_type;
    typedef typename Traits<type>::inverse_type inverse_type;
   
    enum { mcolsize = Traits<type>::msize };
    enum { mrowsize = Traits<type>::msize };
    enum { msize = Traits<type>::msize };
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
    enum { munit = Traits<type>::munit };

    typedef typename Traits<value_type>::real_type real_type;
    typedef typename Traits<value_type>::complex_type complex_type;
    enum { misreal = Traits<value_type>::isreal };
    enum { miscomplex = Traits<value_type>::iscomplex };

    typedef typename Traits<type>::const_row_range_type const_row_range_type;
    typedef typename Traits<type>::const_col_range_type const_col_range_type;
    typedef typename Traits<type>::const_diag_type const_diag_type;
    typedef typename Traits<type>::const_diag_range_type const_diag_range_type;

    typedef typename Traits<type>::const_subtrimatrix_type const_subtrimatrix_type;
    typedef typename Traits<type>::const_subtrimatrix_step_type const_subtrimatrix_step_type;
    typedef typename Traits<type>::const_submatrix_type const_submatrix_type;
    typedef typename Traits<type>::const_submatrix_step_type const_submatrix_step_type;
    typedef typename Traits<type>::const_subvector_type const_subvector_type;

    typedef typename Traits<type>::const_view_type const_view_type;
    typedef typename Traits<type>::const_cview_type const_cview_type;
    typedef typename Traits<type>::const_fview_type const_fview_type;
    typedef typename Traits<type>::const_xview_type const_xview_type;
    typedef typename Traits<type>::const_cmview_type const_cmview_type;
    typedef typename Traits<type>::const_rmview_type const_rmview_type;
    typedef typename Traits<type>::const_transpose_type const_transpose_type;
    typedef typename Traits<type>::const_conjugate_type const_conjugate_type;
    typedef typename Traits<type>::const_adjoint_type const_adjoint_type;

    typedef typename Traits<type>::const_offdiag_type const_offdiag_type;
    typedef typename Traits<type>::const_unitdiag_type const_unitdiag_type;
    typedef typename Traits<type>::const_realview_type const_realview_type;
    typedef typename Traits<type>::const_imagview_type const_imagview_type;
    typedef typename Traits<type>::const_nonconj_type const_nonconj_type;
    typedef typename Traits<type>::nonconst_type nonconst_type;



    //
    // Constructor
    //

    inline BaseMatrix_Tri() {}
    inline BaseMatrix_Tri(const BaseMatrix_Tri<M>&) {}
    inline ~BaseMatrix_Tri() {}

  private:
    void operator=(const BaseMatrix_Tri<M>&);
  public:


    //
    // Access 
    //

    inline const_row_range_type get_row(int i, int j1, int j2) const
    { return const_row_range_type(cptr()+i*stepi()+j1*stepj(),j2-j1,stepj()); }

    inline const_col_range_type get_col(int j, int i1, int i2) const
    { return const_col_range_type(cptr()+j*stepj()+i1*stepi(),i2-i1,stepi()); }

    inline const_diag_range_type get_diag(int i) const
    {
      return const_diag_range_type(
          cptr() + (i<0?(-i*stepi()):(i*stepj())),
          ( i<0 ?  size()+i : size()-i ),
          diagstep());
    }

    inline const_diag_range_type get_diag(int i, int j1, int j2) const
    {
      return const_diag_range_type(
          cptr() + (i<0?(-i*stepi()):(i*stepj())) + j1*diagstep(),
          j2-j1, diagstep());
    }


    // The regular versions respect the indexing style for i and j:
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

    inline const_diag_type diag() const
    {
      CheckDiagTri<mshape>(0);
      return const_diag_type(cptr(),size(),diagstep()); 
    }

    inline const_diag_range_type diag(int i) const
    {
      CheckDiagIndex<mfort>(i,colsize(),rowsize());
      CheckDiagTri<mshape>(i);
      return get_diag(i);
    }

    inline const_diag_range_type diag(int i, int j1, int j2) const
    {
      CheckDiagIndex<mfort>(i,j1,j2,colsize(),rowsize());
      CheckDiagTri<mshape>(i);
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
    inline const_subtrimatrix_type CSubTriMatrix(int i1, int i2) const
    {
      return const_subtrimatrix_type(
          cptr()+i1*diagstep(), i2-i1, stepi(), stepj());
    }

    inline const_subtrimatrix_step_type CSubTriMatrix(
        int i1, int i2, int istep) const
    {
      return const_subtrimatrix_step_type(
          cptr()+i1*diagstep(), (i2-i1)/istep, istep*stepi(), istep*stepj());
    }

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


    // These check the indices according the the indexing style being
    // used, and then calls the above CStyle versions.
    inline const_submatrix_type SubMatrix(
        int i1, int i2, int j1, int j2) const
    {
      CheckRowRange<mfort>(i1,i2,colsize());
      CheckColRange<mfort>(j1,j2,rowsize());
      // Check each corner of submatrix:
      CheckTri<mshape>(i1,j1);
      if (j2 > j1) CheckTri<mshape>(i1,j2-1);
      if (i2 > i1) CheckTri<mshape>(i2-1,j1);
      if (i2 > i1 && j2 > j1) CheckTri<mshape>(i2-1,j2-1);
      return CSubMatrix(i1,i2,j1,j2);
    }

    inline const_submatrix_step_type SubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) const
    {
      CheckRowRange<mfort>(i1,i2,istep,colsize());
      CheckColRange<mfort>(j1,j2,jstep,rowsize());
      // Check each corner of submatrix:
      CheckTri<mshape>(i1,j1);
      if (j2 > j1) CheckTri<mshape>(i1,j2-jstep);
      if (i2 > i1) CheckTri<mshape>(i2-istep,j1);
      if (i2 > i1 && j2 > j1) CheckTri<mshape>(i2-istep,j2-jstep);
      return CSubMatrix(i1,i2,j1,j2,istep,jstep);
    }

    inline const_subvector_type SubVector(
        int i, int j, int istep, int jstep, int s) const
    {
      CheckMatSubVector<mfort>(i,j,istep,jstep,s,colsize(),rowsize());
      // Check each end of subvector:
      CheckTri<mshape>(i,j);
      if (s) CheckTri<mshape>(i+(s-1)*istep,j+(s-1)*jstep);
      return CSubVector(i,j,istep,jstep,s);
    }

    inline const_subtrimatrix_type SubTriMatrix(int i1, int i2) const
    {
      CheckRange<mfort>(i1,i2,size());
      return CSubTriMatrix(i1,i2);
    }

    inline const_subtrimatrix_step_type SubTriMatrix(
        int i1, int i2, int istep) const
    {
      CheckRange<mfort>(i1,i2,istep,size());
      return CSubTriMatrix(i1,i2,istep);
    }


    //
    // Views
    //

    inline const_view_type View() const
    { return const_view_type(cptr(),size(),stepi(),stepj()); }

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
    { return const_transpose_type(cptr(),size(),stepj(),stepi()); }

    inline const_conjugate_type Conjugate() const
    { return const_conjugate_type(cptr(),size(),stepi(),stepj()); }

    inline const_adjoint_type Adjoint() const
    { return const_adjoint_type(cptr(),size(),stepj(),stepi()); }

    inline const_offdiag_type OffDiag(int noff=1) const
    {
      CheckIndex<CStyle>(noff,size());
      return const_offdiag_type(
          cptr()+noff*(this->isupper()?stepj():stepi()),
          size()-noff,stepi(),stepj()); 
    }

    inline const_unitdiag_type ViewAsUnitDiag() const
    { return const_unitdiag_type(cptr(),size(),stepi(),stepj()); }

    inline const_realview_type Real() const
    {
      return const_realview_type(
          reinterpret_cast<const real_type*>(cptr()), size(),
          misreal ? stepi() : 2*stepi(), misreal ? stepj() : 2*stepj());
    }

    inline const_imagview_type Imag() const
    {
      TMVStaticAssert(miscomplex);
      TMVAssert(!munit);
      return const_imagview_type(
          reinterpret_cast<const real_type*>(cptr())+1, size(),
          misreal ? stepi() : 2*stepi(), misreal ? stepj() : 2*stepj());
    }

    inline const_nonconj_type NonConj() const
    { return const_nonconj_type(cptr(),size(),stepi(),stepj()); }

    inline nonconst_type NonConst() const
    { 
      return nonconst_type(
          const_cast<value_type*>(cptr()),size(),stepi(),stepj()); 
    }


    //
    // I/O
    //

    inline void WriteCompact(std::ostream& os) const
    { tmv::WriteCompact(os,CView()); }
    inline void WriteCompact(std::ostream& os, real_type thresh) const
    { tmv::WriteCompact(os,CView(),thresh); }



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
    inline bool isunit() const { return munit; }
    inline bool isupper() const { return ShapeTraits<mshape>::upper; }

    // Note that these last functions need to be defined in a more derived
    // class than this, or an infinite loop will result when compiling.
    // Also, cref, get_row and get_col from BaseMatrix.

    inline size_t colsize() const { return mat().size(); }
    inline size_t rowsize() const { return mat().size(); }
    inline size_t size() const { return mat().size(); }
    inline int stepi() const { return mat().stepi(); }
    inline int stepj() const { return mat().stepj(); }

    inline const value_type* cptr() const { return mat().cptr(); }

  }; // BaseMatrix_Tri

  template <class M>
  class BaseMatrix_Tri_Mutable :
    public BaseMatrix_Tri<M>,
    public BaseMatrix_Mutable<M>
  {
  public:

    typedef M type;
    typedef BaseMatrix_Tri<M> base_tri;

    typedef typename Traits<type>::value_type value_type;
    typedef typename Traits<type>::calc_type calc_type;
    typedef typename Traits<type>::eval_type eval_type;
    typedef typename Traits<type>::copy_type copy_type;

    enum { mcolsize = Traits<type>::msize };
    enum { mrowsize = Traits<type>::msize };
    enum { msize = Traits<type>::msize };
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
    enum { munit = Traits<type>::munit };

    typedef typename Traits<value_type>::real_type real_type;
    typedef typename Traits<value_type>::complex_type complex_type;
    enum { misreal = Traits<value_type>::isreal };
    enum { miscomplex = Traits<value_type>::iscomplex };

    typedef typename Traits<type>::const_row_range_type const_row_range_type;
    typedef typename Traits<type>::const_col_range_type const_col_range_type;
    typedef typename Traits<type>::const_diag_type const_diag_type;
    typedef typename Traits<type>::const_diag_range_type const_diag_range_type;

    typedef typename Traits<type>::const_submatrix_type const_submatrix_type;
    typedef typename Traits<type>::const_submatrix_step_type const_submatrix_step_type;
    typedef typename Traits<type>::const_subvector_type const_subvector_type;
    typedef typename Traits<type>::const_subtrimatrix_type const_subtrimatrix_type;
    typedef typename Traits<type>::const_subtrimatrix_step_type const_subtrimatrix_step_type;

    typedef typename Traits<type>::const_view_type const_view_type;
    typedef typename Traits<type>::const_cview_type const_cview_type;
    typedef typename Traits<type>::const_fview_type const_fview_type;
    typedef typename Traits<type>::const_xview_type const_xview_type;
    typedef typename Traits<type>::const_cmview_type const_cmview_type;
    typedef typename Traits<type>::const_rmview_type const_rmview_type;
    typedef typename Traits<type>::const_transpose_type const_transpose_type;
    typedef typename Traits<type>::const_conjugate_type const_conjugate_type;
    typedef typename Traits<type>::const_adjoint_type const_adjoint_type;

    typedef typename Traits<type>::const_offdiag_type const_offdiag_type;
    typedef typename Traits<type>::const_unitdiag_type const_unitdiag_type;
    typedef typename Traits<type>::const_realview_type const_realview_type;
    typedef typename Traits<type>::const_imagview_type const_imagview_type;
    typedef typename Traits<type>::const_nonconj_type const_nonconj_type;

    typedef typename Traits<type>::row_range_type row_range_type;
    typedef typename Traits<type>::col_range_type col_range_type;
    typedef typename Traits<type>::diag_type diag_type;
    typedef typename Traits<type>::diag_range_type diag_range_type;
    typedef typename Traits<type>::submatrix_type submatrix_type;
    typedef typename Traits<type>::submatrix_step_type submatrix_step_type;
    typedef typename Traits<type>::subvector_type subvector_type;
    typedef typename Traits<type>::subtrimatrix_type subtrimatrix_type;
    typedef typename Traits<type>::subtrimatrix_step_type subtrimatrix_step_type;

    typedef typename Traits<type>::view_type view_type;
    typedef typename Traits<type>::cview_type cview_type;
    typedef typename Traits<type>::fview_type fview_type;
    typedef typename Traits<type>::xview_type xview_type;
    typedef typename Traits<type>::cmview_type cmview_type;
    typedef typename Traits<type>::rmview_type rmview_type;
    typedef typename Traits<type>::transpose_type transpose_type;
    typedef typename Traits<type>::conjugate_type conjugate_type;
    typedef typename Traits<type>::adjoint_type adjoint_type;

    typedef typename Traits<type>::offdiag_type offdiag_type;
    typedef typename Traits<type>::unitdiag_type unitdiag_type;
    typedef typename Traits<type>::realview_type realview_type;
    typedef typename Traits<type>::imagview_type imagview_type;
    typedef typename Traits<type>::nonconj_type nonconj_type;

    typedef typename Traits<type>::reference reference;


    //
    // Constructor
    //

    inline BaseMatrix_Tri_Mutable() {}
    inline BaseMatrix_Tri_Mutable(const BaseMatrix_Tri_Mutable<M>&) {}
    inline ~BaseMatrix_Tri_Mutable() {}


    //
    // Access 
    //

    inline reference operator()(int i, int j)
    {
      CheckRowIndex<mfort>(i,size());
      CheckColIndex<mfort>(j,size());
      CheckTri<mshape>(i,j);
      return ref(i,j);
    }

    inline row_range_type get_row(int i, int j1, int j2) 
    { return row_range_type(ptr()+i*stepi()+j1*stepj(),j2-j1,stepj()); }

    inline col_range_type get_col(int j, int i1, int i2) 
    { return col_range_type(ptr()+j*stepj()+i1*stepi(),i2-i1,stepi()); }

    inline diag_range_type get_diag(int i) 
    {
      return diag_range_type(
          ptr() + (i<0?(-i*stepi()):(i*stepj())),
          ( i<0 ?  size()+i : size()-i ),
          diagstep());
    }

    inline diag_range_type get_diag(int i, int j1, int j2) 
    {
      return diag_range_type(
          ptr() + (i<0?(-i*stepi()):(i*stepj())) + j1*diagstep(),
          j2-j1, diagstep());
    }


    // The regular versions respect the indexing style for i and j:
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

    inline diag_type diag() 
    {
      CheckDiagTri<mshape>(0);
      return diag_type(ptr(),size(),diagstep()); 
    }

    inline diag_range_type diag(int i) 
    {
      CheckDiagIndex<mfort>(i,colsize(),rowsize());
      CheckDiagTri<mshape>(i);
      return get_diag(i);
    }

    inline diag_range_type diag(int i, int j1, int j2) 
    {
      CheckDiagIndex<mfort>(i,j1,j2,colsize(),rowsize());
      CheckDiagTri<mshape>(i);
      return get_diag(i,j1,j2);
    }

    // We need to repeat the const versions so the non-const ones
    // don't clobber them.
    inline value_type operator()(int i, int j) const
    { return base_tri::operator()(i,j); }

    inline const_row_range_type get_row(int i, int j1, int j2) const
    { return base_tri::get_row(i,j1,j2); }
    inline const_col_range_type get_col(int j, int i1, int i2) const
    { return base_tri::get_col(j,i1,i2); }
    inline const_diag_range_type get_diag(int i) const
    { return base_tri::get_diag(i); }
    inline const_diag_range_type get_diag(int i, int j1, int j2) const
    { return base_tri::get_diag(i,j1,j2); }

    inline const_row_range_type row(int i, int j1, int j2) const
    { return base_tri::row(i,j1,j2); }
    inline const_col_range_type col(int j, int i1, int i2) const
    { return base_tri::col(j,i1,i2); }
    inline const_diag_type diag() const
    { return base_tri::diag(); }
    inline const_diag_range_type diag(int i) const
    { return base_tri::diag(i); }
    inline const_diag_range_type diag(int i, int j1, int j2) const
    { return base_tri::diag(i,j1,j2); }


    //
    // Op =
    //

    inline type& operator=(BaseMatrix_Tri_Mutable<M>& m2) 
    {
      TMVStaticAssert((Sizes<msize,M::msize>::same));
      TMVAssert(size() == m2.size());
      TMVStaticAssert((ShapeTraits2<M::mshape,mshape>::assignable));
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

    template <class M2>
    inline type& operator=(const BaseMatrix_Diag<M2>& m2) 
    {
      TMVStaticAssert((Sizes<msize,M2::msize>::same));
      TMVAssert(size() == m2.size());
      m2.diag().AssignTo(diag());
      this->OffDiag().Zero();
      return mat(); 
    }

    inline type& operator=(const value_type x)
    {
      TMVStaticAssert((Sizes<mrowsize,mcolsize>::same));
      TMVAssert(colsize() == rowsize());
      SetToIdentity(x);
      return mat();
    }


    //
    // Modifying Functions
    //

    // These are required for all BaseMatrix_Mutable
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

    // Some more that are added for Tri shape:
    inline type& InvertSelf() 
    { tmv::InvertTri(mat()); return mat(); }

    inline type& SetToIdentity(const value_type x=value_type(1))
    { Zero(); if (!munit) diag().SetAllTo(x); return mat(); }


    //
    // SubMatrix, etc.
    //

    // These versions always uses CStyle
    inline subtrimatrix_type CSubTriMatrix(int i1, int i2) 
    {
      return subtrimatrix_type(ptr()+i1*diagstep(), i2-i1, stepi(), stepj()); 
    }

    inline subtrimatrix_step_type CSubTriMatrix(int i1, int i2, int istep) 
    {
      return subtrimatrix_step_type(
          ptr()+i1*diagstep(), (i2-i1)/istep, istep*stepi(), istep*stepj());
    }

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

    // These check the indices according the the indexing style being
    // used, and then calls the above CStyle versions.
    inline submatrix_type SubMatrix(int i1, int i2, int j1, int j2) 
    {
      CheckRowRange<mfort>(i1,i2,colsize());
      CheckColRange<mfort>(j1,j2,rowsize());
      // Check each corner of submatrix:
      CheckTri<mshape>(i1,j1);
      if (j2 > j1) CheckTri<mshape>(i1,j2-1);
      if (i2 > i1) CheckTri<mshape>(i2-1,j1);
      if (i2 > i1 && j2 > j1) CheckTri<mshape>(i2-1,j2-1);
      return CSubMatrix(i1,i2,j1,j2);
    }

    inline submatrix_step_type SubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) 
    {
      CheckRowRange<mfort>(i1,i2,istep,colsize());
      CheckColRange<mfort>(j1,j2,jstep,rowsize());
      // Check each corner of submatrix:
      CheckTri<mshape>(i1,j1);
      if (j2 > j1) CheckTri<mshape>(i1,j2-jstep);
      if (i2 > i1) CheckTri<mshape>(i2-istep,j1);
      if (i2 > i1 && j2 > j1) CheckTri<mshape>(i2-istep,j2-jstep);
      return CSubMatrix(i1,i2,j1,j2,istep,jstep);
    }

    inline subvector_type SubVector(
        int i, int j, int istep, int jstep, int s) 
    {
      CheckMatSubVector<mfort>(i,j,istep,jstep,s,colsize(),rowsize());
      // Check each end of subvector:
      CheckTri<mshape>(i,j);
      if (s) CheckTri<mshape>(i+(s-1)*istep,j+(s-1)*jstep);
      return CSubVector(i,j,istep,jstep,s);
    }

    inline subtrimatrix_type SubTriMatrix(int i1, int i2) 
    {
      CheckRange<mfort>(i1,i2,size());
      return CSubTriMatrix(i1,i2);
    }

    inline subtrimatrix_step_type SubTriMatrix(
        int i1, int i2, int istep) 
    {
      CheckRange<mfort>(i1,i2,istep,size());
      return CSubTriMatrix(i1,i2,istep);
    }


    // Repeat the const versions:
    inline const_submatrix_type CSubMatrix(
        int i1, int i2, int j1, int j2) const
    { return base_tri::CSubMatrix(i1,i2,j1,j2); }
    inline const_submatrix_step_type CSubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) const
    { return base_tri::CSubMatrix(i1,i2,j1,j2,istep,jstep); }
    inline const_subvector_type CSubVector(
        int i, int j, int istep, int jstep, int s) const
    { return base_tri::CSubVector(i,j,istep,jstep,s); }
    inline const_subtrimatrix_type CSubTriMatrix(int i1, int i2) const
    { return base_tri::CSubTriMatrix(i1,i2); }
    inline const_subtrimatrix_step_type CSubTriMatrix(
        int i1, int i2, int istep) const
    { return base_tri::CSubTriMatrix(i1,i2,istep); }

    inline const_submatrix_type SubMatrix(
        int i1, int i2, int j1, int j2) const
    { return base_tri::SubMatrix(i1,i2,j1,j2); }
    inline const_submatrix_step_type SubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) const
    { return base_tri::SubMatrix(i1,i2,j1,j2,istep,jstep); }
    inline const_subvector_type SubVector(
        int i, int j, int istep, int jstep, int s) const
    { return base_tri::SubVector(i,j,istep,jstep,s); }
    inline const_subtrimatrix_type SubTriMatrix(int i1, int i2) const
    { return base_tri::SubTriMatrix(i1,i2); }
    inline const_subtrimatrix_step_type SubTriMatrix(
        int i1, int i2, int istep) const
    { return base_tri::SubTriMatrix(i1,i2,istep); }


    //
    // Views
    //

    inline view_type View() 
    { return view_type(ptr(),size(),stepi(),stepj()); }

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
    { return transpose_type(ptr(),size(),stepj(),stepi()); }

    inline conjugate_type Conjugate() 
    { return conjugate_type(ptr(),size(),stepi(),stepj()); }

    inline adjoint_type Adjoint() 
    { return adjoint_type(ptr(),size(),stepj(),stepi()); }

    inline offdiag_type OffDiag(int noff=1) 
    {
      CheckIndex<CStyle>(noff,size());
      return offdiag_type(
          ptr()+noff*(this->isupper()?stepj():stepi()),
          size()-noff,stepi(),stepj()); 
    }

    inline unitdiag_type ViewAsUnitDiag()
    { return const_unitdiag_type(ptr(),size(),stepi(),stepj()); }

    inline realview_type Real() 
    {
      return realview_type(reinterpret_cast<real_type*>(ptr()), size(),
          misreal ? stepi() : 2*stepi(), misreal ? stepj() : 2*stepj());
    }

    inline imagview_type Imag() 
    {
      TMVStaticAssert(miscomplex);
      TMVAssert(!munit);
      return imagview_type(reinterpret_cast<real_type*>(ptr())+1, size(),
          misreal ? stepi() : 2*stepi(), misreal ? stepj() : 2*stepj());
    }

    inline nonconj_type NonConj()
    { return nonconj_type(ptr(),size(),stepi(),stepj()); }


    // Repeat the const versions:
    inline const_view_type View() const
    { return base_tri::View(); }
    inline const_cview_type CView() const
    { return base_tri::CView(); }
    inline const_fview_type FView() const
    { return base_tri::FView(); }
    inline const_xview_type XView() const
    { return base_tri::XView(); }
    inline const_cmview_type CMView() const
    { return base_tri::XView(); }
    inline const_rmview_type RMView() const
    { return base_tri::XView(); }
    inline const_transpose_type Transpose() const
    { return base_tri::Transpose(); }
    inline const_conjugate_type Conjugate() const
    { return base_tri::Conjugate(); }
    inline const_adjoint_type Adjoint() const
    { return base_tri::Adjoint(); }
    inline const_offdiag_type OffDiag(int noff=1) const
    { return base_tri::OffDiag(noff); }
    inline const_unitdiag_type ViewAsUnitDiag() const
    { return base_tri::ViewAsUnitDiag(); }
    inline const_realview_type Real() const
    { return base_tri::Real(); }
    inline const_imagview_type Imag() const
    { return base_tri::Imag(); }
    inline const_nonconj_type NonConj() const
    { return base_tri::NonConj(); }


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

    inline int diagstep() const
    { return mdiagstep == UNKNOWN ? stepi() + stepj() : mdiagstep; }

    // Note that these last functions need to be defined in a more derived
    // class than this, or an infinite loop will result when compiling.
    // Also, cref and cptr from above.

    inline size_t colsize() const { return mat().size(); }
    inline size_t rowsize() const { return mat().size(); }
    inline size_t size() const { return mat().size(); }
    inline int stepi() const { return mat().stepi(); }
    inline int stepj() const { return mat().stepj(); }

    inline value_type* ptr() { return mat().ptr(); }
    inline reference ref(int i, int j) { return mat().ref(i,j); }
    inline value_type cref(int i, int j) { return mat().cref(i,j); }

  }; // BaseMatrix_Tri_Mutable

  //
  // Zero
  //

  template <int algo, class M> struct ZeroU_Helper;

  // algo 1: LowerTri -- Transpose
  template <class M> 
  struct ZeroU_Helper<1,M> 
  {
    static inline void call(M& m) 
    { m.Transpose().Zero(); } 
  };

  // algo 2: RowMajor NonUnit
  template <class M> 
  struct ZeroU_Helper<2,M> 
  {
    static void call(M& m) 
    {
      const int n = m.size();
      for(int i=0;i<n;++i) m.get_row(i,i,n).Zero();
    } 
  };

  // algo 3: ColMajor NonUnit
  template <class M> 
  struct ZeroU_Helper<3,M> 
  {
    static void call(M& m) 
    {
      const int n = m.size();
      for(int j=0;j<n;++j) m.get_col(j,0,j+1).Zero();
    } 
  };

  // algo 4: RowMajor Unit
  template <class M> 
  struct ZeroU_Helper<4,M> 
  {
    static void call(M& m) 
    {
      const int n = m.size();
      for(int i=0;i<n;++i) m.get_row(i,i+1,n).Zero();
    } 
  };

  // algo 5: ColMajor Unit
  template <class M> 
  struct ZeroU_Helper<5,M> 
  {
    static void call(M& m) 
    {
      const int n = m.size();
      for(int j=0;j<n;++j) m.get_col(j,0,j).Zero();
    } 
  };

  template <class M>
  inline void Zero(BaseMatrix_Tri_Mutable<M>& m)
  {
    enum { algo = (
        Shape(M::mshape) == LowerTri ? 1 :
        M::munit ? M::mrowmajor ? 4 : 5 :
        M::mrowmajor ? 2 : 3 ) };
    ZeroU_Helper<algo,M>::call(m.mat());
  }

  //
  // SetAllTo
  //

  template <int algo, class M, class T> struct SetAllToU_Helper;

  // algo 1: LowerTri -- Transpose
  template <class M, class T> 
  struct SetAllToU_Helper<1,M,T>
  { 
    static inline void call(M& m, const T& val) 
    { m.Transpose().SetAllTo(val); } 
  };

  // algo 2: RowMajor NonUnit
  template <class M, class T> 
  struct SetAllToU_Helper<2,M,T>
  {
    static void call(M& m, const T& val) 
    {
      const int n = m.size();
      for(int i=0;i<n;++i) m.get_row(i,i,n).SetAllTo(val);
    } 
  };

  // algo 3: ColMajor
  template <class M, class T> 
  struct SetAllToU_Helper<3,M,T> 
  {
    static void call(M& m, const T& val) 
    {
      const int n = m.size();
      for(int j=0;j<n;++j) m.get_col(j,0,j+1).SetAllTo(val);
    } 
  };

  // algo 4: RowMajor Unit
  template <class M, class T> 
  struct SetAllToU_Helper<4,M,T> 
  {
    static void call(M& m, const T& val) 
    {
      const int n = m.size();
      for(int i=0;i<n;++i) m.get_row(i,i+1,n).SetAllTo(val);
    } 
  };

  // algo 5: ColMajor Unit
  template <class M, class T> 
  struct SetAllToU_Helper<5,M,T> 
  {
    static void call(M& m, const T& val) 
    {
      const int n = m.size();
      for(int j=0;j<n;++j) m.get_col(j,0,j).SetAllTo(val);
    } 
  };

  template <class M, class T>
  inline void SetAllTo(BaseMatrix_Tri_Mutable<M>& m, const T& val)
  {
    enum { algo = (
        Shape(M::mshape) == LowerTri ? 1 :
        M::munit ? M::mrowmajor ? 4 : 5 :
        M::mrowmajor ? 2 : 3 ) };
    SetAllToU_Helper<algo,M,T>::call(m.mat(),val);
  }

  //
  // AddToAll
  //

  template <int algo, class M, class T> struct AddToAllU_Helper;

  // algo 1: LowerTri -- Transpose
  template <class M, class T> 
  struct AddToAllU_Helper<1,M,T> 
  { 
    static inline void call(M& m, const T& val) 
    { m.Transpose().AddToAll(val); } 
  };

  // algo 2: RowMajor NonUnit
  template <class M, class T> 
  struct AddToAllU_Helper<2,M,T> 
  {
    static void call(M& m, const T& val) 
    {
      const int n = m.size();
      for(int i=0;i<n;++i) m.get_row(i,i,n).AddToAll(val);
    } 
  };

  // algo 3: ColMajor
  template <class M, class T> 
  struct AddToAllU_Helper<3,M,T> 
  {
    static void call(M& m, const T& val) 
    {
      const int n = m.size();
      for(int j=0;j<n;++j) m.get_col(j,0,j+1).AddToAll(val);
    } 
  };

  // algo 4: RowMajor Unit
  template <class M, class T> 
  struct AddToAllU_Helper<4,M,T> 
  {
    static void call(M& m, const T& val) 
    {
      const int n = m.size();
      for(int i=0;i<n;++i) m.get_row(i,i+1,n).AddToAll(val);
    } 
  };

  // algo 5: ColMajor Unit
  template <class M, class T> 
  struct AddToAllU_Helper<5,M,T> 
  {
    static void call(M& m, const T& val) 
    {
      const int n = m.size();
      for(int j=0;j<n;++j) m.get_col(j,0,j).AddToAll(val);
    } 
  };

  template <class M, class T>
  inline void AddToAll(BaseMatrix_Tri_Mutable<M>& m, const T& val)
  {
    enum { algo = (
        Shape(M::mshape) == LowerTri ? 1 :
        M::munit ? M::mrowmajor ? 4 : 5 :
        M::mrowmajor ? 2 : 3 ) };
    AddToAllU_Helper<algo,M,T>::call(m.mat(),val);
  }

  //
  // Clip
  //

  template <int algo, class M, class RT> struct ClipU_Helper;

  // algo 1: LowerTri -- Transpose
  template <class M, class RT> 
  struct ClipU_Helper<1,M,RT>
  {
    static inline void call(M& m, const RT& thresh) 
    { m.Transpose().Clip(thresh); } 
  };

  // algo 2: RowMajor NonUnit
  template <class M, class RT> 
  struct ClipU_Helper<2,M,RT> 
  {
    static void call(M& m, const RT& thresh) 
    {
      const int n = m.size();
      for(int i=0;i<n;++i) m.get_row(i,i,n).Clip(thresh);
    } 
  };

  // algo 3: ColMajor
  template <class M, class RT> 
  struct ClipU_Helper<3,M,RT> 
  {
    static void call(M& m, const RT& thresh) 
    {
      const int n = m.size();
      for(int j=0;j<n;++j) m.get_col(j,0,j+1).Clip(thresh);
    } 
  };

  // algo 4: RowMajor Unit
  template <class M, class RT> 
  struct ClipU_Helper<4,M,RT> 
  {
    static void call(M& m, const RT& thresh) 
    {
      const int n = m.size();
      for(int i=0;i<n;++i) m.get_row(i,i+1,n).Clip(thresh);
    } 
  };

  // algo 5: ColMajor Unit
  template <class M, class RT> 
  struct ClipU_Helper<5,M,RT> 
  {
    static void call(M& m, const RT& thresh) 
    {
      const int n = m.size();
      for(int j=0;j<n;++j) m.get_col(j,0,j).Clip();
    } 
  };

  template <class M, class RT>
  inline void Clip(BaseMatrix_Tri_Mutable<M>& m, const RT& thresh)
  {
    enum { algo = (
        Shape(M::mshape) == LowerTri ? 1 :
        M::munit ? M::mrowmajor ? 4 : 5 :
        M::mrowmajor ? 2 : 3 ) };
    ClipU_Helper<algo,M,RT>::call(m.mat(),thresh);
  }

  //
  // ApplyToAll
  //

  template <int algo, class M, class F> struct ApplyToAllU_Helper;

  // algo 1: LowerTri -- Transpose
  template <class M, class F> 
  struct ApplyToAllU_Helper<1,M,F>
  { 
    static inline void call(M& m, const F& f) 
    { m.Transpose().ApplyToAll(f); } 
  };

  // algo 2: RowMajor NonUnit
  template <class M, class F> 
  struct ApplyToAllU_Helper<2,M,F> 
  {
    static void call(M& m, const F& f) 
    {
      const int n = m.size();
      for(int i=0;i<n;++i) m.get_row(i,i,n).ApplyToAll(f);
    } 
  };

  // algo 3: ColMajor
  template <class M, class F> 
  struct ApplyToAllU_Helper<3,M,F> 
  {
    static void call(M& m, const F& f) 
    {
      const int n = m.size();
      for(int j=0;j<n;++j) m.get_col(j,0,j+1).ApplyToAll(f);
    } 
  };

  // algo 4: RowMajor Unit
  template <class M, class F> 
  struct ApplyToAllU_Helper<4,M,F> 
  {
    static void call(M& m, const F& f) 
    {
      const int n = m.size();
      for(int i=0;i<n;++i) m.get_row(i,i+1,n).ApplyToAll(f);
    } 
  };

  // algo 5: ColMajor Unit
  template <class M, class F> 
  struct ApplyToAllU_Helper<5,M,F> 
  {
    static void call(M& m, const F& f) 
    {
      const int n = m.size();
      for(int j=0;j<n;++j) m.get_col(j,0,j).ApplyToAll();
    } 
  };

  template <class M, class F>
  inline void ApplyToAll(BaseMatrix_Tri_Mutable<M>& m, const F& f)
  {
    enum { algo = (
        Shape(M::mshape) == LowerTri ? 1 :
        M::munit ? M::mrowmajor ? 4 : 5 :
        M::mrowmajor ? 2 : 3 ) };
    ApplyToAllU_Helper<algo,M,F>::call(m.mat(),f);
  }

  //
  // ConjugateSelf
  //

  template <int algo, class M> struct ConjugateU_Helper;

  // algo 0: Not complex, nothing to do
  template <class M>
  struct ConjugateU_Helper<0,M>
  { static inline void call(M& ) {} };

  // algo 1: LowerTri -- Transpose
  template <class M>
  struct ConjugateU_Helper<1,M>
  {
    static inline void call(M& m)
    { m.Transpose().ConjugateSelf(); }
  };

  // algo 2: RowMajor NonUnit
  template <class M>
  struct ConjugateU_Helper<2,M> 
  {
    static void call(M& m)
    {
      const int n = m.size();
      for(int i=0;i<n;++i) m.get_row(i,i,n).ConjugateSelf();
    } 
  };

  // algo 3: ColMajor
  template <class M>
  struct ConjugateU_Helper<3,M> 
  {
    static void call(M& m)
    {
      const int n = m.size();
      for(int j=0;j<n;++j) m.get_col(j,0,j+1).ConjugateSelf();
    } 
  };

  // algo 4: RowMajor Unit
  template <class M> 
  struct ConjugateU_Helper<4,M> 
  {
    static void call(M& m)
    {
      const int n = m.size();
      for(int i=0;i<n;++i) m.get_row(i,i+1,n).ConjugateSelf();
    } 
  };

  // algo 5: ColMajor Unit
  template <class M> 
  struct ConjugateU_Helper<5,M> 
  {
    static void call(M& m)
    {
      const int n = m.size();
      for(int j=0;j<n;++j) m.get_col(j,0,j).ConjugateSelf();
    } 
  };

  template <class M>
  inline void ConjugateSelf(BaseMatrix_Tri_Mutable<M>& m)
  {
    enum { algo = (
        Shape(M::mshape) == LowerTri ? 1 :
        M::munit ? M::mrowmajor ? 4 : 5 :
        M::mrowmajor ? 2 : 3 ) };
    ConjugateU_Helper<algo,M>::call(m.mat());
  }

  //
  // Copy Matrices
  //

  // Defined in TMV_CopyU.h
  template <class M1, class M2> 
  inline void Copy(
      const BaseMatrix_Tri<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2);

  template <bool upper, class M1, class M2>
  struct CopyUM_Helper;
  template <class M1, class M2>
  struct CopyUM_Helper<true,M1,M2> // m1 is uppertri
  {
    static void call(const M1& m1, M2& m2)
    {
      m2.UpperTri() = m1;
      m2.LowerTri().OffDiag().Zero();
    }
  };
  template <class M1, class M2>
  struct CopyUM_Helper<false,M1,M2> // m1 is lowertri
  {
    static void call(const M1& m1, M2& m2)
    {
      m2.LowerTri() = m1;
      m2.UpperTri().OffDiag().Zero();
    }
  };
  template <class M1, class M2> 
  inline void Copy(
      const BaseMatrix_Tri<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
  {
    enum { upper = ShapeTraits<M1::mshape>::upper };
    CopyUM_Helper<upper,M1,M2>::call(m1.mat(),m2.mat()); 
  }

  //
  // Arithmetic that trivially translates into simpler operations:
  //
  // M = x * U
  //

  template <bool upper, int ix, class T, class M1, class M2>
  struct MultXUM_Helper;
  template <int ix, class T, class M1, class M2>
  struct MultXUM_Helper<true,ix,T,M1,M2> // m1 is uppertri
  {
    static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    {
      m2.UpperTri() = x*m1;
      m2.LowerTri().OffDiag().Zero();
    }
  };
  template <int ix, class T, class M1, class M2>
  struct MultXUM_Helper<false,ix,T,M1,M2> // m1 is lowertri
  {
    static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    {
      m2.LowerTri() = x*m1;
      m2.UpperTri().OffDiag().Zero();
    }
  };
  template <int ix, class T, class M1, class M2>
  inline void MultXM(
      const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
      BaseMatrix_Rec_Mutable<M2>& m2)
  {
    enum { upper = ShapeTraits<M1::mshape>::upper };
    MultXUM_Helper<upper,ix,T,M1,M2>::call(x,m1.mat(),m2.mat());
  }


  //
  // M += x * U
  //

  template <bool upper, int ix, class T, class M1, class M2>
  struct AddUM_Helper;
  template <int ix, class T, class M1, class M2>
  struct AddUM_Helper<true,ix,T,M1,M2> // m1 is uppertri
  {
    static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    { m2.UpperTri() += x*m1; }
  };
  template <int ix, class T, class M1, class M2>
  struct AddUM_Helper<false,ix,T,M1,M2> // m1 is lowertri
  {
    static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
    { m2.LowerTri() += x*m1; }
  };
  template <int ix, class T, class M1, class M2>
  inline void AddMM(
      const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
      BaseMatrix_Rec_Mutable<M2>& m2)
  {
    enum { upper = ShapeTraits<M1::mshape>::upper };
    AddUM_Helper<upper,ix,T,M1,M2>::call(x,m1.mat(),m2.mat());
  }


  //
  //
  // U += x * D
  //

  template <int ix, class T, class M1, class M2>
  static void AddMM(
      const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
      BaseMatrix_Tri_Mutable<M2>& m2)
  { m2.diag() += x * m1.diag(); }



  //
  // TypeText 
  //

  template <class M>
  static std::string TypeText(const BaseMatrix_Tri<M>& m)
  {
    std::ostringstream s;
    s << "BaseMatrix_Tri< "<<TypeText(m.mat())<<" >";
    return s.str();
  }

  template <class M>
  static std::string TypeText(const BaseMatrix_Tri_Mutable<M>& m)
  {
    std::ostringstream s;
    s << "BaseMatrix_Tri_Mutable< "<<TypeText(m.mat())<<" >";
    return s.str();
  }

} // namespace tmv

#endif
