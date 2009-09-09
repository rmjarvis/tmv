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
// This file defines the BaseMatrix_Diag and BaseMatrix_Diag_Mutable classes.
//
// See TMV_DiagMatrix.h for the functions that are defined for these objects.
//

#ifndef TMV_BaseMatrix_Diag_H
#define TMV_BaseMatrix_Diag_H

#include "tmv/TMV_BaseMatrix.h"
#include "tmv/TMV_MultVV.h"
#include "tmv/TMV_BaseMatrix_Rec.h"

namespace tmv {

  template <class M>
  class BaseMatrix_Diag;
  template <class M>
  class BaseMatrix_Diag_Mutable;

  template <class M>
  static void Read(std::istream& is, BaseMatrix_Diag_Mutable<M>& m);

  template <class M1, class M2> 
  static void Copy(
      const BaseMatrix_Diag<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2);
  template <class M1, class M2> 
  static void Copy(
      const BaseMatrix_Diag<M1>& m1, BaseMatrix_Diag_Mutable<M2>& m2);

  template <class M>
  class BaseMatrix_Diag :
    public BaseMatrix_Calc<M>
  {
  public:
    typedef M type;

    typedef typename Traits<type>::value_type value_type;
    typedef typename Traits<value_type>::real_type real_type;
    typedef typename Traits<value_type>::complex_type complex_type;
    enum { misreal = Traits<value_type>::isreal };
    enum { miscomplex = Traits<value_type>::iscomplex };

    typedef typename Traits<type>::calc_type calc_type;
    typedef typename Traits<type>::eval_type eval_type;
    typedef typename Traits<type>::copy_type copy_type;

    enum { mcolsize = Traits<type>::msize };
    enum { mrowsize = Traits<type>::msize };
    enum { msize = Traits<type>::msize };
    enum { mfort = Traits<type>::mfort };
    enum { mshape = Traits<type>::mshape };
    enum { mrowmajor = Traits<type>::mrowmajor };
    enum { mcolmajor = Traits<type>::mcolmajor };
    enum { mcalc = true };
    enum { mstep = Traits<type>::mstep };
    enum { mconj = Traits<type>::mconj };

    typedef typename Traits<type>::const_diag_type const_diag_type;

    typedef typename Traits<type>::const_subdiagmatrix_type const_subdiagmatrix_type;
    typedef typename Traits<type>::const_subdiagmatrix_step_type const_subdiagmatrix_step_type;
    typedef typename Traits<type>::const_view_type const_view_type;
    typedef typename Traits<type>::const_cview_type const_cview_type;
    typedef typename Traits<type>::const_fview_type const_fview_type;
    typedef typename Traits<type>::const_xview_type const_xview_type;
    typedef typename Traits<type>::const_cmview_type const_cmview_type;
    typedef typename Traits<type>::const_rmview_type const_rmview_type;
    typedef typename Traits<type>::const_transpose_type const_transpose_type;
    typedef typename Traits<type>::const_conjugate_type const_conjugate_type;
    typedef typename Traits<type>::const_adjoint_type const_adjoint_type;
    typedef typename Traits<type>::const_realview_type const_realview_type;
    typedef typename Traits<type>::const_imagview_type const_imagview_type;
    typedef typename Traits<type>::const_nonconj_type const_nonconj_type;
    typedef typename Traits<type>::nonconst_type nonconst_type;



    //
    // Constructor
    //

    inline BaseMatrix_Diag() {}
    inline BaseMatrix_Diag(const BaseMatrix_Diag<M>&) {}
    inline ~BaseMatrix_Diag() {}

  private:
    void operator=(const BaseMatrix_Diag<M>&);
  public:


    //
    // Access 
    //

    inline value_type operator()(int i, int j) const
    {
      CheckRowIndex<mfort>(i,size());
      CheckColIndex<mfort>(j,size());
      TMVAssert(i == j);
      return cref(i);
    }

    inline value_type operator()(int i) const
    {
      CheckIndex<mfort>(i,size());
      return cref(i);
    }

    inline const_diag_type diag() const
    { return const_diag_type(cptr(),size(),step()); }


    //
    // Functions
    //

    inline value_type SumElements() const
    { return diag().SumElements(); }

    inline real_type SumAbsElements() const
    { return diag().SumAbsElements(); }

    inline real_type MaxAbsElement() const
    { return diag().MaxAbsElement(); }

    inline real_type NormSq() const
    { return diag().NormSq(); }

    inline real_type NormSq(const real_type scale) const
    { return diag().NormSq(scale); }

    inline real_type NormF() const 
    { return diag().Norm2(); }

    inline real_type Norm() const
    { return NormF(); }

    inline real_type Norm1() const
    { return diag().MaxAbsElement(); }

    inline real_type Norm2() const
    { return diag().MaxAbsElement(); }

    inline real_type NormInf() const
    { return diag().MaxAbsElement(); }

    template <class ret_type, class F>
    inline ret_type SumElements(const F& f) const
    { return diag().SumElements(f); }



    //
    // SubMatrix, etc.
    //

    // These versions always uses CStyle
    inline const_subdiagmatrix_type CSubDiagMatrix(int i1, int i2) const
    {
      return const_subdiagmatrix_type(
          cptr()+i1*step(), i2-i1, step());
    }

    inline const_subdiagmatrix_step_type CSubDiagMatrix(
        int i1, int i2, int istep) const
    {
      return const_subdiagmatrix_step_type(
          cptr()+i1*step(), (i2-i1)/istep, istep*step());
    }


    // These check the indices according the the indexing style being
    // used, and then calls the above CStyle versions.
    inline const_subdiagmatrix_type SubDiagMatrix(int i1, int i2) const
    {
      CheckRange<mfort>(i1,i2,size());
      return CSubDiagMatrix(i1,i2);
    }

    inline const_subdiagmatrix_step_type SubDiagMatrix(
        int i1, int i2, int istep) const
    {
      CheckRange<mfort>(i1,i2,istep,size());
      return CSubDiagMatrix(i1,i2,istep);
    }


    //
    // Views
    //

    inline const_view_type View() const
    { return const_view_type(cptr(),size(),step()); }

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
    { return View(); }

    inline const_conjugate_type Conjugate() const
    { return const_conjugate_type(cptr(),size(),step()); }

    inline const_adjoint_type Adjoint() const
    { return Conjugate(); }

    inline const_realview_type Real() const
    {
      return const_realview_type(
          reinterpret_cast<const real_type*>(cptr()), size(),
          misreal ? step() : 2*step());
    }

    inline const_imagview_type Imag() const
    {
      TMVStaticAssert(miscomplex);
      return const_imagview_type(
          reinterpret_cast<const real_type*>(cptr())+1, size(),
          misreal ? step() : 2*step());
    }

    inline const_nonconj_type NonConj() const
    { return const_nonconj_type(cptr(),size(),step()); }

    inline nonconst_type NonConst() const
    { return nonconst_type(const_cast<value_type*>(cptr()),size(),step()); }


    //
    // I/O
    //

    inline void WriteCompact(std::ostream& os) const
    { os << "D " << diag() << std::endl; }
    inline void WriteCompact(std::ostream& os, real_type thresh) const
    { os << "D "; diag().Write(os,thresh); os << std::endl; }


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

    inline bool isconj() const { return mconj; }

    // Note that these last functions need to be defined in a more derived
    // class than this, or an infinite loop will result.
    // Also, cref, get_row and get_col from BaseMatrix.

    inline size_t colsize() const { return mat().size(); }
    inline size_t rowsize() const { return mat().size(); }
    inline size_t size() const { return mat().size(); }
    inline int step() const { return mat().step(); }

    inline const value_type* cptr() const { return mat().cptr(); }
    inline value_type cref(int i, int j) const 
    { return (i!=j ? value_type(0) : cref(i)); }
    inline value_type cref(int i) const 
    { return mat().cref(i); }

  };

  // Specify ExactSameStorage for diagonal matrices:
  template <class M>
  inline bool ExactSameStorage(
      const BaseMatrix_Diag<M>& m1, const BaseMatrix_Diag<M>& m2)
  { return m1.step() == m2.step(); }

  template <class M>
  class BaseMatrix_Diag_Mutable :
    public BaseMatrix_Diag<M>,
    public BaseMatrix_Mutable<M>
  {
  public:

    typedef M type;
    typedef BaseMatrix_Diag<M> base_diag;

    typedef typename Traits<type>::value_type value_type;
    typedef typename Traits<value_type>::real_type real_type;
    typedef typename Traits<value_type>::complex_type complex_type;
    enum { misreal = Traits<value_type>::isreal };
    enum { miscomplex = Traits<value_type>::iscomplex };

    typedef typename Traits<type>::calc_type calc_type;
    typedef typename Traits<type>::eval_type eval_type;
    typedef typename Traits<type>::copy_type copy_type;

    enum { mcolsize = Traits<type>::msize };
    enum { mrowsize = Traits<type>::msize };
    enum { msize = Traits<type>::msize };
    enum { mfort = Traits<type>::mfort };
    enum { mshape = Traits<type>::mshape };
    enum { mrowmajor = Traits<type>::mrowmajor };
    enum { mcolmajor = Traits<type>::mcolmajor };
    enum { mcalc = true };
    enum { mstep = Traits<type>::mstep };
    enum { mconj = Traits<type>::mconj };

    typedef typename Traits<type>::const_diag_type const_diag_type;

    typedef typename Traits<type>::const_subdiagmatrix_type const_subdiagmatrix_type;
    typedef typename Traits<type>::const_subdiagmatrix_step_type const_subdiagmatrix_step_type;
    typedef typename Traits<type>::const_view_type const_view_type;
    typedef typename Traits<type>::const_cview_type const_cview_type;
    typedef typename Traits<type>::const_fview_type const_fview_type;
    typedef typename Traits<type>::const_xview_type const_xview_type;
    typedef typename Traits<type>::const_cmview_type const_cmview_type;
    typedef typename Traits<type>::const_rmview_type const_rmview_type;
    typedef typename Traits<type>::const_transpose_type const_transpose_type;
    typedef typename Traits<type>::const_conjugate_type const_conjugate_type;
    typedef typename Traits<type>::const_adjoint_type const_adjoint_type;
    typedef typename Traits<type>::const_realview_type const_realview_type;
    typedef typename Traits<type>::const_imagview_type const_imagview_type;
    typedef typename Traits<type>::const_nonconj_type const_nonconj_type;

    typedef typename Traits<type>::diag_type diag_type;
    typedef typename Traits<type>::subdiagmatrix_type subdiagmatrix_type;
    typedef typename Traits<type>::subdiagmatrix_step_type subdiagmatrix_step_type;
    typedef typename Traits<type>::view_type view_type;
    typedef typename Traits<type>::cview_type cview_type;
    typedef typename Traits<type>::fview_type fview_type;
    typedef typename Traits<type>::xview_type xview_type;
    typedef typename Traits<type>::cmview_type cmview_type;
    typedef typename Traits<type>::rmview_type rmview_type;
    typedef typename Traits<type>::transpose_type transpose_type;
    typedef typename Traits<type>::conjugate_type conjugate_type;
    typedef typename Traits<type>::adjoint_type adjoint_type;
    typedef typename Traits<type>::realview_type realview_type;
    typedef typename Traits<type>::imagview_type imagview_type;
    typedef typename Traits<type>::nonconj_type nonconj_type;

    typedef typename Traits<type>::reference reference;


    //
    // Constructor
    //

    inline BaseMatrix_Diag_Mutable() {}
    inline BaseMatrix_Diag_Mutable(const BaseMatrix_Diag_Mutable<M>&) {}
    inline ~BaseMatrix_Diag_Mutable() {}


    //
    // Access 
    //

    inline reference operator()(int i, int j)
    {
      TMVAssert(i == j);
      CheckIndex<mfort>(i,size());
      return ref(i);
    }

    inline reference operator()(int i)
    {
      CheckIndex<mfort>(i,size());
      return ref(i);
    }

    inline diag_type diag() 
    { return diag_type(ptr(),size(),step()); }


    // We need to repeat the const versions so the non-const ones
    // don't clobber them.
    inline value_type operator()(int i, int j) const
    { return base_diag::operator()(i,j); }
    inline const_diag_type diag() const
    { return base_diag::diag(); }


    //
    // Op =
    //

    inline type& operator=(BaseMatrix_Diag_Mutable<M>& m2) 
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


    //
    // Modifying Functions
    //

    inline type& Zero() 
    { return SetAllTo(value_type(0)); }

    inline type& SetAllTo(value_type x) 
    { diag().SetAllTo(x); return mat(); }

    inline type& AddToAll(value_type x) 
    { diag().AddToAll(x); return mat(); }

    inline type& Clip(real_type thresh) 
    { diag().Clip(thresh); return mat(); }

    template <class F>
    inline type& ApplyToAll(const F& f)
    { diag().ApplyToAll(f); return mat(); }

    inline type& ConjugateSelf() 
    { diag().ConjugateSelf(); return mat(); }

    inline type& TransposeSelf() 
    { return mat(); }

    inline type& InvertSelf() 
    {
      const int n=size();
      for(int i=0;i<n;++i) ref(i) = real_type(1)/cref(i);
      return mat(); 
    }

    inline type& SetToIdentity(const value_type x=value_type(1))
    { diag().SetAllTo(x); return mat(); }

    inline type& Swap(int i1, int i2) 
    {
      CheckIndex<mfort>(i1,size());
      CheckIndex<mfort>(i2,size());
      diag().CSwap(i1,i2);
      return mat();
    }

    inline type& CPermute(const int*const p, int i1, int i2) 
    { diag().Permute(p,i1,i2); return mat(); }
    inline type& Permute(const int*const p, int i1, int i2) 
    {
      CheckRange<mfort>(i1,i2,size());
      return CPermute(p,i1,i2);
    }
    inline type& Permute(const int*const p) 
    { CPermute(p,0,size()); return mat(); }

    inline type& CReversePermute(const int*const p, int i1, int i2) 
    { diag().ReversePermute(p,i1,i2); }
    inline type& ReversePermute(const int*const p, int i1, int i2) 
    {
      CheckRange<mfort>(i1,i2,size());
      return CReversePermute(p,i1,i2);
    }
    inline type& ReversePermute(const int*const p) 
    { CReversePermute(p,0,size()); return mat(); }


    //
    // SubDiagMatrix, etc.
    //

    // These versions always uses CStyle
    inline subdiagmatrix_type CSubDiagMatrix(int i1, int i2) 
    { return subdiagmatrix_type(ptr()+i1*step(), i2-i1, step()); }

    inline subdiagmatrix_step_type CSubDiagMatrix(int i1, int i2, int istep) 
    {
      return subdiagmatrix_step_type(
          ptr()+i1*step(), (i2-i1)/istep, istep*step());
    }

    // These check the indices according the the indexing style being
    // used, and then calls the above CStyle versions.
    inline subdiagmatrix_type SubDiagMatrix(int i1, int i2) 
    {
      CheckRange<mfort>(i1,i2,size());
      return CSubDiagMatrix(i1,i2);
    }

    inline subdiagmatrix_step_type SubMatrix(
        int i1, int i2, int istep) 
    {
      CheckRange<mfort>(i1,i2,istep,size());
      return CSubDiagMatrix(i1,i2,istep);
    }


    // Repeat the const versions:
    inline const_subdiagmatrix_type CSubDiagMatrix(int i1, int i2) const
    { return base_diag::CSubDiagMatrix(i1,i2); }
    inline const_subdiagmatrix_step_type CSubDiagMatrix(
        int i1, int i2, int istep) const
    { return base_diag::CSubDiagMatrix(i1,i2,istep); }

    inline const_subdiagmatrix_type SubDiagMatrix(int i1, int i2) const
    { return base_diag::SubDiagMatrix(i1,i2); }
    inline const_subdiagmatrix_step_type SubDiagMatrix(
        int i1, int i2, int istep) const
    { return base_diag::SubDiagMatrix(i1,i2,istep); }


    //
    // Views
    //

    inline view_type View() 
    { return view_type(ptr(),size(),step()); }

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
    { return View(); }

    inline conjugate_type Conjugate() 
    { return conjugate_type(ptr(),size(),step()); }

    inline adjoint_type Adjoint() 
    { return Conjugate(); }

    inline realview_type Real() 
    {
      return realview_type(reinterpret_cast<real_type*>(ptr()), size(),
          misreal ? step() : 2*step());
    }

    inline imagview_type Imag() 
    {
      TMVStaticAssert(miscomplex);
      return imagview_type(reinterpret_cast<real_type*>(ptr())+1, size(),
          misreal ? step() : 2*step());
    }

    inline nonconj_type NonConj()
    { return nonconj_type(ptr(),size(),step()); }


    // Repeat the const versions:
    inline const_view_type View() const
    { return base_diag::View(); }
    inline const_cview_type CView() const
    { return base_diag::CView(); }
    inline const_fview_type FView() const
    { return base_diag::FView(); }
    inline const_xview_type XView() const
    { return base_diag::XView(); }
    inline const_cmview_type CMView() const
    { return base_diag::XView(); }
    inline const_rmview_type RMView() const
    { return base_diag::XView(); }
    inline const_transpose_type Transpose() const
    { return base_diag::Transpose(); }
    inline const_conjugate_type Conjugate() const
    { return base_diag::Conjugate(); }
    inline const_adjoint_type Adjoint() const
    { return base_diag::Adjoint(); }
    inline const_realview_type Real() const
    { return base_diag::Real(); }
    inline const_imagview_type Imag() const
    { return base_diag::Imag(); }
    inline const_nonconj_type NonConj() const
    { return base_diag::NonConj(); }



    // Defined in TMV_DiagMatrixIO.h
    inline void Read(std::istream& is)
    { tmv::Read(is,mat()); }


    //
    // Auxilliary routines
    //

    inline const type& mat() const
    { return *static_cast<const type*>(this); }
    inline type& mat()
    { return *static_cast<type*>(this); }


    // Note that these last functions need to be defined in a more derived
    // class than this, or an infinite loop will result when compiling.
    // Also, cref and cptr from above.

    inline size_t colsize() const { return mat().size(); }
    inline size_t rowsize() const { return mat().size(); }
    inline size_t size() const { return mat().size(); }
    inline int step() const { return mat().step(); }

    inline value_type* ptr() { return mat().ptr(); }
    inline reference ref(int i) { return mat().ref(i); }
    inline value_type cref(int i) { return mat().cref(i); }

  }; // BaseMatrix_Diag_Mutable

  template <class T, IndexStyle I=CStyle>
  class DiagMatrix;
  template <class T, int S=UNKNOWN, bool C=false, IndexStyle I=CStyle>
  class ConstDiagMatrixView;
  template <class T, int S=UNKNOWN, bool C=false, IndexStyle I=CStyle>
  class DiagMatrixView;
  template <class T, int N, IndexStyle I=CStyle>
  class SmallDiagMatrix;
  template <class T, int N, int S, bool C=false, IndexStyle I=CStyle>
  class ConstSmallDiagMatrixView;
  template <class T, int N, int S, bool C=false, IndexStyle I=CStyle>
  class SmallDiagMatrixView;

  // This helper class helps decide calc_type for composite classes:
  template <class T, int cs, int rs, bool rm, bool fort>
  struct MCopyHelper<T,Diag,cs,rs,rm,fort>
  { typedef SmallDiagMatrix<T,cs,fort?FortranStyle:CStyle> type; };
  template <class T, int rs, bool rm, bool fort>
  struct MCopyHelper<T,Diag,UNKNOWN,rs,rm,fort>
  { typedef SmallDiagMatrix<T,rs,fort?FortranStyle:CStyle> type; };
  template <class T, int cs, bool rm, bool fort>
  struct MCopyHelper<T,Diag,cs,UNKNOWN,rm,fort>
  { typedef SmallDiagMatrix<T,cs,fort?FortranStyle:CStyle> type; };
  template <class T, bool rm, bool fort>
  struct MCopyHelper<T,Diag,UNKNOWN,UNKNOWN,rm,fort>
  { typedef DiagMatrix<T,fort?FortranStyle:CStyle> type; };

  //
  // Copy Matrices
  //

  template <class M1, class M2> 
  static void Copy(
      const BaseMatrix_Diag<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
  {
    if (SameStorage(m1,m2)) {
      m2.diag() = m1.diag();
      m2.UpperTri().OffDiag().Zero();
      m2.LowerTri().OffDiag().Zero();
    }
    else {
      m2.Zero();
      m2.diag() = m1.diag();
    }
  }

  template <class M1, class M2> 
  inline void Copy(
      const BaseMatrix_Diag<M1>& m1, BaseMatrix_Diag_Mutable<M2>& m2)
  { m2.diag() = m1.diag(); }


  //
  // Swap Matrices
  //

  template <class M1, class M2> 
  inline void Swap(
      BaseMatrix_Diag_Mutable<M1>& m1, BaseMatrix_Diag_Mutable<M2>& m2)
  { Swap(m1.diag(),m2.diag()); }


  //
  // Matrix ==, != Matrix
  //

  template <class M1, class M2>
  inline bool operator==(
      const BaseMatrix_Diag<M1>& m1, const BaseMatrix_Diag<M2>& m2)
  {
    TMVStaticAssert((Sizes<M1::mcolsize,M2::mcolsize>::same)); 
    TMVStaticAssert((Sizes<M1::mrowsize,M2::mrowsize>::same)); 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    return m1.diag() == m2.diag();
  }

  template <class M1, class M2>
  inline bool operator!=(
      const BaseMatrix_Diag<M1>& m1, const BaseMatrix_Diag<M2>& m2)
  { return !(m1 == m2); }


  // Most of the diagmatrix arithmetic trivially translates into a
  // vector arithmetic calculation.
  // The only exception is Matrix * DiagMatrix, so that one is in MultMD.h.  
  // The rest are here:

  //
  // D += x * D
  // M += x * D
  //

  template <int ix1, class T1, class M1, class M2>
  inline void AddMM(
      const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
      BaseMatrix_Diag_Mutable<M2>& m2)
  { m2.diag() += x1 * m1.diag(); }

  template <int ix1, class T1, class M1, class M2>
  inline void AddMM(
      const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
      BaseMatrix_Rec_Mutable<M2>& m2)
  { m2.diag() += x1 * m1.diag(); }

  //
  // D = x * D + x * D
  // M = x * D + x * D
  // M = x * D + x * M
  // M = x * M + x * D
  //

  template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  inline void AddMM(
      const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
      const Scaling<ix2,T2>& x2, const BaseMatrix_Diag<M2>& m2,
      BaseMatrix_Diag_Mutable<M3>& m3)
  { m3.diag() = x1 * m1.diag() + x2 * m2.diag(); }

  template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  inline void InlineAddMM(
      const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
      const Scaling<ix2,T2>& x2, const BaseMatrix_Diag<M2>& m2,
      BaseMatrix_Rec_Mutable<M3>& m3)
  { m3.Zero(); m3.diag() = x1 * m1.diag() + x2 * m2.diag(); }

  template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  inline void InlineAddMM(
      const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
      const Scaling<ix2,T2>& x2, const BaseMatrix_Calc<M2>& m2,
      BaseMatrix_Rec_Mutable<M3>& m3)
  { (m3 = x2 * m2).diag() += x1 * m1.diag(); }

  template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
  inline void InlineAddMM(
      const Scaling<ix1,T1>& x1, const BaseMatrix_Calc<M1>& m1,
      const Scaling<ix2,T2>& x2, const BaseMatrix_Diag<M2>& m2,
      BaseMatrix_Rec_Mutable<M3>& m3)
  { (m3 = x1 * m1).diag() += x2 * m2.diag(); }

  //
  // AddDX
  //

  template <class T, class M>
  inline void AddMX(const T& x2, BaseMatrix_Diag_Mutable<M>& m3)
  { m3.diag().AddToAll(x2); }

  template <int ix1, class T1, class M1, class T2, class M3>
  inline void AddMX(const Scaling<ix1,T1>& x1, const BaseMatrix<M1>& m1,
      const T2& x2, BaseMatrix_Diag_Mutable<M3>& m3)
  {
    TMVStaticAssert((Sizes<M1::mcolsize,M1::mrowsize>::same));
    TMVAssert(m1.colsize() == m1.rowsize());
    AddMX(x2,(m3 = x1 * m1));
  }


  //
  // MultXD
  //

  template <int ix, class T, class M>
  inline void MultXM(
      const Scaling<ix,T>& x, BaseMatrix_Diag_Mutable<M>& m)
  { m.diag() *= x; }

  template <int ix1, class T1, class M1, class M2>
  inline void MultXM(
      const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
      BaseMatrix_Diag_Mutable<M2>& m2)
  { m2.diag() = x1 * m1.diag(); }

  template <int ix1, class T1, class M1, class M2>
  static void MultXM(
      const Scaling<ix1,T1>& x1, const BaseMatrix_Diag<M1>& m1,
      BaseMatrix_Rec_Mutable<M2>& m2)
  {
    if (SameStorage(m1,m2)) {
      m2.diag() = x1 * m1.diag();
      m2.UpperTri().OffDiag().Zero();
      m2.LowerTri().OffDiag().Zero();
    }
    else {
      m2.Zero();
      m2.diag() = x1 * m1.diag();
    }
  }


  // 
  // MultDV
  //

  template <bool add, int ix, class T, class M1, class V2, class V3>
  inline void MultMV(const Scaling<ix,T>& x,
      const BaseMatrix_Diag<M1>& m1, const BaseVector_Calc<V2>& v2,
      BaseVector_Mutable<V3>& v3)
  { ElemMultVV<add>(x,m1.diag(),v2.vec(),v3.vec()); }

  template <bool add, int ix, class T, class V1, class M2, class V3>
  inline void MultVM(const Scaling<ix,T>& x,
      const BaseVector_Calc<V1>& v1, const BaseMatrix_Diag<M2>& m2,
      BaseVector_Mutable<V3>& v3)
  { ElemMultVV<add>(x,m2.diag(),v1.vec(),v3.vec()); }

  template <class V1, int ix, class T, class M2>
  inline void MultEqVM(BaseVector_Mutable<V1>& v1,
      const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
  { ElemMultVV<false>(x,m2.diag(),v1.vec(),v1.vec()); }


  //
  // MultDD
  //

  template <bool add, int ix, class T, class M1, class M2, class M3>
  inline void MultMM(const Scaling<ix,T>& x,
      const BaseMatrix_Diag<M1>& m1, const BaseMatrix_Diag<M2>& m2,
      BaseMatrix_Diag_Mutable<M3>& m3)
  { 
    typename M3::diag_type m3d = m3.diag();
    ElemMultVV<add>(x,m1.diag(),m2.diag(),m3d); 
  }

  template <class M1, int ix, class T, class M2>
  inline void MultEqMM(BaseMatrix_Diag_Mutable<M1>& m1,
      const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
  { 
    typename M1::diag_type m1d = m1.diag();
    ElemMultVV<false>(x,m1.diag(),m2.diag(),m1d); 
  }

  template <bool add, int ix, class T, class M1, class M2, class M3>
  static void MultMM(const Scaling<ix,T>& x,
      const BaseMatrix_Diag<M1>& m1, const BaseMatrix_Diag<M2>& m2,
      BaseMatrix_Rec_Mutable<M3>& m3)
  {
    typename M3::diag_type m3d = m3.diag();
    if (!add && SameStorage(m1,m3) || SameStorage(m2,m3)) {
      ElemMultVV<add>(x,m1.diag(),m2.diag(),m3d); 
      m3.UpperTri().OffDiag().Zero();
      m3.LowerTri().OffDiag().Zero();
    }
    else {
      if (!add) m3.Zero();
      ElemMultVV<add>(x,m1.diag(),m2.diag(),m3d); 
    }
  }


  //
  // TypeText 
  //

  template <class M>
  static std::string TypeText(const BaseMatrix_Diag<M>& m)
  {
    std::ostringstream s;
    s << "BaseMatrix_Diag< "<<TypeText(m.mat())<<" >";
    return s.str();
  }

  template <class M>
  static std::string TypeText(const BaseMatrix_Diag_Mutable<M>& m)
  {
    std::ostringstream s;
    s << "BaseMatrix_Diag_Mutable< "<<TypeText(m.mat())<<" >";
    return s.str();
  }

} // namespace tmv

#endif
