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
// This file defines the TMV SmallTriMatrix class.
//
// Constructors:
//
//    As with the comments in TMV_TriMatrix.h, I omit the Upper or Lower
//    before TriMatrix in the below constructors, since they are the same
//    form for each.
//
//    SmallTriMatrix<T,dt,stor,I>(size_t n)
//        Makes a Triangular Matrix with column size = row size = n
//        with _uninitialized_ values.
//
//    SmallTriMatrix<T,dt,stor,I>(size_t n, T x)
//        Makes a Triangular Matrix with column size = row size = n
//        with all values = x
//
//    SmallTriMatrix<T,dt,stor,I>(size_t n, T* vv)
//    SmallTriMatrix<T,dt,stor,I>(size_t n, const std::vector<T>& vv)
//        Makes a Triangular Matrix with column size = row size = n
//        which copies the values from vv.
//
//    SmallTriMatrix<T,dt,stor,I>(const Matrix<T>& m)
//    SmallTriMatrix<T,dt,stor,I>(const SmallTriMatrix<T>& m)
//        Makes a TriMatrix which copies the corresponding elements of m.
//
//

#ifndef TMV_SmallTriMatrix_H
#define TMV_SmallTriMatrix_H

#include "tmv/TMV_BaseMatrix_Tri.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_SmallMatrix.h"
#include <vector>

namespace tmv {

  template <class T, int N, DiagType D, StorageType S, IndexStyle I> 
  struct Traits<SmallUpperTriMatrix<T,N,D,S,I> >
  {
    typedef T value_type;

    typedef typename Traits<T>::real_type real_type;
    typedef typename Traits<T>::complex_type complex_type;
    enum { misreal = Traits<T>::isreal };
    enum { miscomplex = Traits<T>::iscomplex };

    typedef SmallUpperTriMatrix<T,N,D,S,I> type;
    typedef const type& calc_type;
    typedef const type& eval_type;
    typedef type copy_type;

    enum { mcolsize = N };
    enum { mrowsize = N };
    enum { msize = N };
    enum { mshape = UpperTri };
    enum { mfort = (I == FortranStyle) };
    enum { mcalc = true };
    enum { mrowmajor = (S == RowMajor) };
    enum { mcolmajor = (S == ColMajor) };
    enum { mstor = S };
    enum { mstepi = (S==RowMajor ? N : 1) };
    enum { mstepj = (S==RowMajor ? 1 : N) };
    enum { mdiagstep = N+1 };
    enum { mconj = false };
    enum { munit = (D == UnitDiag) };

    enum { twoSi = misreal ? int(mstepi) : int(IntTraits<mstepi>::twoS) };
    enum { twoSj = misreal ? int(mstepj) : int(IntTraits<mstepj>::twoS) };
    enum { notC = miscomplex };

    typedef ConstVectorView<T,mstepi,false,I> const_col_range_type;
    typedef ConstVectorView<T,mstepj,false,I> const_row_range_type;
    typedef ConstSmallVectorView<T,N,mdiagstep,false,I> const_diag_type;
    typedef ConstVectorView<T,mdiagstep,false,I> const_diag_range_type;

    typedef ConstUpperTriMatrixView<T,D,mstepi,mstepj,false,I> const_subtrimatrix_type;
    typedef ConstUpperTriMatrixView<T,D,UNKNOWN,UNKNOWN,false,I> const_subtrimatrix_step_type;
    typedef ConstMatrixView<T,mstepi,mstepj,false,I> const_submatrix_type;
    typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,false,I> const_submatrix_step_type;
    typedef ConstVectorView<T,UNKNOWN,false,I> const_subvector_type;

    typedef ConstSmallUpperTriMatrixView<T,N,D,mstepi,mstepj,false,I> const_view_type;
    typedef ConstSmallUpperTriMatrixView<T,N,D,mstepi,mstepj,false,CStyle> const_cview_type;
    typedef ConstSmallUpperTriMatrixView<T,N,D,mstepi,mstepj,false,FortranStyle> const_fview_type;
    typedef ConstUpperTriMatrixView<T,D> const_xview_type;
    typedef ConstSmallUpperTriMatrixView<T,N,D,1,mstepj,false,I> const_cmview_type;
    typedef ConstSmallUpperTriMatrixView<T,N,D,mstepi,1,false,I> const_rmview_type;
    typedef ConstSmallUpperTriMatrixView<T,N,D,mstepi,mstepj,notC,I> const_conjugate_type;
    typedef ConstSmallLowerTriMatrixView<T,N,D,mstepj,mstepi,false,I> const_transpose_type;
    typedef ConstSmallLowerTriMatrixView<T,N,D,mstepj,mstepi,notC,I> const_adjoint_type;

    typedef ConstUpperTriMatrixView<T,NonUnitDiag,mstepi,mstepj,false,I> const_offdiag_type;
    typedef ConstSmallUpperTriMatrixView<T,N,UnitDiag,mstepi,mstepj,false,I> const_unitdiag_type;
    typedef ConstSmallUpperTriMatrixView<real_type,N,D,twoSi,twoSj,false,I> const_realview_type;
    typedef const_realview_type const_imagview_type;
    typedef ConstSmallUpperTriMatrixView<T,N,D,mstepi,mstepj,false,I> const_nonconj_type;

    typedef InvalidType inverse_type;
    typedef SmallUpperTriMatrixView<T,N,D,mstepi,mstepj,false,I> nonconst_type;

    typedef T& reference;

    typedef VectorView<T,mstepi,false,I> col_range_type;
    typedef VectorView<T,mstepj,false,I> row_range_type;
    typedef VectorView<T,mdiagstep,false,I> diag_type;
    typedef VectorView<T,mdiagstep,false,I> diag_range_type;

    typedef UpperTriMatrixView<T,D,mstepi,mstepj,false,I> subtrimatrix_type;
    typedef UpperTriMatrixView<T,D,UNKNOWN,UNKNOWN,false,I> subtrimatrix_step_type;
    typedef MatrixView<T,mstepi,mstepj,false,I> submatrix_type;
    typedef MatrixView<T,UNKNOWN,UNKNOWN,false,I> submatrix_step_type;
    typedef VectorView<T,UNKNOWN,false,I> subvector_type;

    typedef SmallUpperTriMatrixView<T,N,D,mstepi,mstepj,false,I> view_type;
    typedef SmallUpperTriMatrixView<T,N,D,mstepi,mstepj,false,CStyle> cview_type;
    typedef SmallUpperTriMatrixView<T,N,D,mstepi,mstepj,false,FortranStyle> fview_type;
    typedef UpperTriMatrixView<T,D> xview_type;
    typedef SmallUpperTriMatrixView<T,N,D,1,mstepj,false,I> cmview_type;
    typedef SmallUpperTriMatrixView<T,N,D,mstepi,1,false,I> rmview_type;
    typedef SmallUpperTriMatrixView<T,N,D,mstepi,mstepj,notC,I> conjugate_type;
    typedef SmallLowerTriMatrixView<T,N,D,mstepj,mstepi,false,I> transpose_type;
    typedef SmallLowerTriMatrixView<T,N,D,mstepj,mstepi,notC,I> adjoint_type;

    typedef UpperTriMatrixView<T,NonUnitDiag,mstepi,mstepj,false,I> offdiag_type;
    typedef SmallUpperTriMatrixView<T,N,UnitDiag,mstepi,mstepj,false,I> unitdiag_type;
    typedef SmallUpperTriMatrixView<real_type,N,D,twoSi,twoSj,false,I> realview_type;
    typedef realview_type imagview_type;
    typedef SmallUpperTriMatrixView<T,N,D,mstepi,mstepj,false,I> nonconj_type;
  };

#ifdef XTEST
#ifdef TMVDEBUG
#define XTEST_DEBUG
#endif
#endif

  template <class T, int N, DiagType D, StorageType S, IndexStyle I> 
  class SmallUpperTriMatrix : 
    public BaseMatrix_Tri_Mutable<SmallUpperTriMatrix<T,N,D,S,I> >
  {
  public:

    typedef SmallUpperTriMatrix<T,N,D,S,I> type;
    typedef BaseMatrix_Tri_Mutable<type> base_mut;
    enum { misreal = Traits<type>::misreal };
    enum { miscomplex = Traits<type>::miscomplex };
    enum { mshape = Traits<type>::mshape };
    enum { munit = Traits<type>::munit };
    enum { mstepi = Traits<type>::mstepi };
    enum { mstepj = Traits<type>::mstepj };

    //
    // Constructors
    //

    inline SmallUpperTriMatrix()
    {
      TMVStaticAssert(N>0);
      TMVStaticAssert(S==RowMajor || S==ColMajor); 
#ifdef TMVDEBUG
      this->SetAllTo(T(888));
#endif
    }

    explicit inline SmallUpperTriMatrix(T x)
    {
      TMVStaticAssert(N>0);
      TMVStaticAssert(S==RowMajor || S==ColMajor);
      this->SetAllTo(x);
    }

    explicit inline SmallUpperTriMatrix(const T* vv) 
    {
      TMVStaticAssert(N>0);
      TMVStaticAssert(S==RowMajor || S==ColMajor);
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      ConstSmallMatrixView<T,N,D,mstepi,mstepj> m2(vv);
      SmallMatrixView<T,N,D,mstepi,mstepj> m1(ptr());
      m2.AssignTo(m1);
    }

    explicit inline SmallUpperTriMatrix(const std::vector<T>& vv) 
    {
      TMVStaticAssert(N>0);
      TMVStaticAssert(S==RowMajor || S==ColMajor);
      TMVAssert(vv.size() == N*N);
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      std::copy(vv.begin(),vv.end(),itsm.get());
    }

    template <class M2> 
    inline SmallUpperTriMatrix(const BaseMatrix<M2>& m2) 
    { 
      TMVStaticAssert(N>0);
      TMVStaticAssert(S==RowMajor || S==ColMajor); 
      TMVStaticAssert(M2::misreal || miscomplex);
      TMVStaticAssert((ShapeTraits2<M2::mshape,mshape>::assignable));
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      m2.AssignTo(*this);
    }

    inline SmallUpperTriMatrix(const type& m2) 
    { 
      TMVStaticAssert(N>0);
      TMVStaticAssert(S==RowMajor || S==ColMajor); 
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      m2.AssignTo(*this);
    }

    template <class M2>
    inline SmallUpperTriMatrix(const BaseMatrix_Tri<M2>& m2) 
    {
      TMVStaticAssert(N>0);
      TMVStaticAssert(S==RowMajor || S==ColMajor);
      TMVStaticAssert(Shape(M2::mshape) == UpperTri);
      TMVStaticAssert(M2::munit || !munit);
      TMVStaticAssert(M2::misreal || miscomplex);
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      m2.AssignTo(*this);
      if (isunit() && !m2.isunit()) {
        if (m2.size() > 0) this->OffDiag() = m2.OffDiag();
      }
      else m2.AssignTo(*this);
    }

    template <class M2>
    inline SmallUpperTriMatrix(const BaseMatrix_Rec<M2>& m2) 
    {
      TMVStaticAssert(N>0);
      TMVStaticAssert(S==RowMajor || S==ColMajor);
      TMVStaticAssert(M2::misreal || miscomplex);
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      if (isunit()) {
        m2.UnitUpperTri().AssignTo(*this);
      } else {
        m2.UpperTri().AssignTo(*this);
      }
    }

    inline ~SmallUpperTriMatrix()
    {
#ifdef TMVDEBUG
      this->SetAllTo(T(999));
#endif
    }


    //
    // Op=
    //

    inline type& operator=(const type& m2)
    { 
      if (&m2 != this) 
        base_mut::operator=(m2);
      return *this;
    }

    template <class M2>
    inline type& operator=(const BaseMatrix<M2>& m2)
    {
      base_mut::operator=(m2);
      return *this;
    }

    template <class M2>
    inline type& operator=(const BaseMatrix_Diag<M2>& m2)
    {
      base_mut::operator=(m2);
      return *this;
    }

    inline type& operator=(T x)
    {
      base_mut::operator=(x);
      return *this;
    }


    // 
    // Auxilliary Functions
    //

    inline const T* cptr() const { return itsm; }
    inline T* ptr() { return itsm; }

    inline T cref(int i, int j) const 
    {
      return (
          (D==UnitDiag && i==j ) ? T(1) :
          (i>j) ? T(0) :
          itsm[S==RowMajor ? i*stepi() + j : i + j*stepj()]);
    }

    inline T& ref(int i, int j)
    { return itsm[S==RowMajor ? i*stepi() + j : i + j*stepj()]; }

    inline size_t size() const { return N; }
    inline int stepi() const { return mstepi; }
    inline int stepj() const { return mstepj; }
    inline DiagType dt() const { return D; }
    inline bool isunit() const { return D == UnitDiag; }
    inline bool isconj() const { return false; }
    inline bool isrm() const { return S==RowMajor; }
    inline bool iscm() const { return S==ColMajor; }
    inline StorageType stor() const { return S; }

  protected :

    StackArray<T,N*N> itsm;

  }; // SmallUpperTriMatrix

  template <class T, int N, DiagType D, StorageType S>
  class SmallUpperTriMatrixF : 
    public SmallUpperTriMatrix<T,N,D,S,FortranStyle>
  {
  public:

    typedef SmallUpperTriMatrixF<T,N,D,S> type;
    typedef SmallUpperTriMatrix<T,N,D,S,FortranStyle> mtype;

    inline SmallUpperTriMatrixF() : mtype() {}
    explicit inline SmallUpperTriMatrixF(T x) : mtype(x) {}
    explicit inline SmallUpperTriMatrixF(const T* vv) : mtype(vv) {}
    explicit inline SmallUpperTriMatrixF(const std::vector<T>& vv) : 
      mtype(vv) {}
    template <class M2> 
    inline SmallUpperTriMatrixF(const BaseMatrix<M2>& m2) : mtype(m2) {}
    inline SmallUpperTriMatrixF(const type& m2) : mtype(m2) {}
    template <class M2>
    inline SmallUpperTriMatrixF(const BaseMatrix_Tri<M2>& m2) : mtype(m2) {}
    template <class M2>
    inline SmallUpperTriMatrixF(const BaseMatrix_Rec<M2>& m2) : mtype(m2) {}
    inline ~SmallUpperTriMatrixF() {}

    inline type& operator=(const type& m2)
    { mtype::operator=(m2); return *this; }
    template <class M2>
    inline type& operator=(const BaseMatrix<M2>& m2)
    { mtype::operator=(m2); return *this; }
    template <class M2>
    inline type& operator=(const BaseMatrix_Diag<M2>& m2)
    { mtype::operator=(m2); return *this; }
    inline type& operator=(T x)
    { mtype::operator=(x); return *this; }

  }; // SmallUpperTriMatrixF

  template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I>
  struct Traits<ConstSmallUpperTriMatrixView<T,N,D,Si,Sj,C,I> >
  {
    typedef T value_type;

    typedef typename Traits<T>::real_type real_type;
    typedef typename Traits<T>::complex_type complex_type;
    enum { misreal = Traits<T>::isreal };
    enum { miscomplex = Traits<T>::iscomplex };

    typedef ConstSmallUpperTriMatrixView<T,N,D,Si,Sj,C,I> type;
    typedef const type& calc_type;
    typedef const type& eval_type;
    typedef InvalidType inverse_type;

    enum { mcolsize = N };
    enum { mrowsize = N };
    enum { msize = N };
    enum { mshape = UpperTri };
    enum { mfort = (I == FortranStyle) };
    enum { mcalc = true };
    enum { mrowmajor = (Sj == 1) };
    enum { mcolmajor = (Si == 1) };
    enum { mstor = (mrowmajor ? RowMajor : ColMajor) };
    enum { mstepi = Si };
    enum { mstepj = Sj };
    enum { mdiagstep = IntTraits2<Si,Sj>::sum };
    enum { mconj = C };
    enum { munit = (D == UnitDiag) };

    // In case N == UNKNOWN
    typedef typename MCopyHelper<T,UpperTri,N,N,mrowmajor,mfort>::type copy_type;

    enum { twoSi = misreal ? Si : IntTraits<Si>::twoS };
    enum { twoSj = misreal ? Sj : IntTraits<Sj>::twoS };
    enum { notC = !C && miscomplex };

    typedef ConstVectorView<T,mstepi,C,I> const_col_range_type;
    typedef ConstVectorView<T,mstepj,C,I> const_row_range_type;
    typedef ConstVectorView<T,mdiagstep,C,I> const_diag_type;
    typedef ConstVectorView<T,mdiagstep,C,I> const_diag_range_type;

    typedef ConstUpperTriMatrixView<T,D,mstepi,mstepj,C,I> const_subtrimatrix_type;
    typedef ConstUpperTriMatrixView<T,D,UNKNOWN,UNKNOWN,C,I> const_subtrimatrix_step_type;
    typedef ConstMatrixView<T,mstepi,mstepj,C,I> const_submatrix_type;
    typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,C,I> const_submatrix_step_type;
    typedef ConstVectorView<T,UNKNOWN,C,I> const_subvector_type;

    typedef ConstSmallUpperTriMatrixView<T,N,D,mstepi,mstepj,C,I> const_view_type;
    typedef ConstSmallUpperTriMatrixView<T,N,D,mstepi,mstepj,C,CStyle> const_cview_type;
    typedef ConstSmallUpperTriMatrixView<T,N,D,mstepi,mstepj,C,FortranStyle> const_fview_type;
    typedef ConstUpperTriMatrixView<T,D,UNKNOWN,UNKNOWN,C> const_xview_type;
    typedef ConstSmallUpperTriMatrixView<T,N,D,1,mstepj,C,I> const_cmview_type;
    typedef ConstSmallUpperTriMatrixView<T,N,D,mstepi,1,C,I> const_rmview_type;
    typedef ConstSmallUpperTriMatrixView<T,N,D,mstepi,mstepj,notC,I> const_conjugate_type;
    typedef ConstSmallLowerTriMatrixView<T,N,D,mstepj,mstepi,C,I> const_transpose_type;
    typedef ConstSmallLowerTriMatrixView<T,N,D,mstepj,mstepi,notC,I> const_adjoint_type;

    typedef ConstUpperTriMatrixView<T,NonUnitDiag,mstepi,mstepj,C,I> const_offdiag_type;
    typedef ConstSmallUpperTriMatrixView<T,N,UnitDiag,mstepi,mstepj,C,I> const_unitdiag_type;
    typedef ConstSmallUpperTriMatrixView<real_type,N,D,twoSi,twoSj,false,I> const_realview_type;
    typedef const_realview_type const_imagview_type;
    typedef ConstSmallUpperTriMatrixView<T,N,D,mstepi,mstepj,false,I> const_nonconj_type;
    typedef SmallUpperTriMatrixView<T,N,D,mstepi,mstepj,C,I> nonconst_type;
  };

  template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I>
  class ConstSmallUpperTriMatrixView :
    public BaseMatrix_Tri<ConstSmallUpperTriMatrixView<T,N,D,Si,Sj,C,I> >
  {
  public:

    typedef ConstSmallUpperTriMatrixView<T,N,D,Si,Sj,C,I> type;
    enum { mrowmajor = Traits<type>::mrowmajor };
    enum { mcolmajor = Traits<type>::mcolmajor };
    enum { mshape = Traits<type>::mshape };
    enum { munit = Traits<type>::munit };

    //
    // Constructors
    //
    inline ConstSmallUpperTriMatrixView(const T* m, size_t s, int si, int sj) :
      itsm(m), itss(s), itssi(si), itssj(sj) {}

    inline ConstSmallUpperTriMatrixView(const T* m, size_t s, int si) :
      itsm(m), itss(s), itssi(si), itssj(Sj)
    { TMVStaticAssert(Sj != UNKNOWN); }

    inline ConstSmallUpperTriMatrixView(const T* m, size_t s) :
      itsm(m), itss(s), itssi(Si), itssj(Sj)
    { TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); }

    inline ConstSmallUpperTriMatrixView(const T* m) :
      itsm(m), itss(N), itssi(Si), itssj(Sj)
    {
      TMVStaticAssert(N != UNKNOWN);
      TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
    }

    inline ConstSmallUpperTriMatrixView(const type& m2) :
      itsm(m2.cptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int Si2, int Sj2, IndexStyle I2>
    inline ConstSmallUpperTriMatrixView(
        const ConstUpperTriMatrixView<T,D,Si2,Sj2,C,I2>& m2) :
      itsm(m2.cptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int Si2, int Sj2, IndexStyle I2>
    inline ConstSmallUpperTriMatrixView(
        const UpperTriMatrixView<T,D,Si2,Sj2,C,I2>& m2) :
      itsm(m2.cptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int N2, int Si2, int Sj2, IndexStyle I2>
    inline ConstSmallUpperTriMatrixView(
        const ConstSmallUpperTriMatrixView<T,N2,D,Si2,Sj2,C,I2>& m2) :
      itsm(m2.cptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int N2, int Si2, int Sj2, IndexStyle I2>
    inline ConstSmallUpperTriMatrixView(
        const SmallUpperTriMatrixView<T,N2,D,Si2,Sj2,C,I2>& m2) :
      itsm(m2.cptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    inline ~ConstSmallUpperTriMatrixView() {
#ifdef TMV_DEBUG
      itsm = 0;
#endif
    }

  private :
    inline void operator=(const type& m2);
  public :

    //
    // Auxilliary Functions
    //

    inline const T* cptr() const { return itsm; }

    inline T cref(int i, int j) const 
    {
      return (
          (D==UnitDiag && i==j ) ? T(1) :
          (i>j) ? T(0) :
          DoConj<C>(itsm[i*stepi() + j*stepj()]));
    }

    inline size_t colsize() const { return itss; }
    inline size_t rowsize() const { return itss; }
    inline size_t size() const { return itss; }
    inline int stepi() const { return itssi; }
    inline int stepj() const { return itssj; }
    inline bool isconj() const { return C; }
    inline bool isunit() const { return D == UnitDiag; }
    inline bool isrm() const
    { return mrowmajor || (!mcolmajor && stepj() == 1); }
    inline bool iscm() const
    { return mcolmajor || (!mrowmajor && stepi() == 1); }

  private :

#ifdef TMV_DEBUG
    const T* itsm;
#else
    const T*const itsm;
#endif
    const StepInt<N> itss;
    const StepInt<Si> itssi;
    const StepInt<Sj> itssj;

  }; // ConstSmallUpperTriMatrixView

  template <class T, int N, DiagType D, int Si, int Sj, bool C>
  class ConstSmallUpperTriMatrixViewF :
    public ConstSmallUpperTriMatrixView<T,N,D,Si,Sj,C,FortranStyle>
  {
  public:

    typedef ConstSmallUpperTriMatrixViewF<T,N,D,Si,Sj,C> type;
    typedef ConstSmallUpperTriMatrixView<T,N,D,Si,Sj,C,FortranStyle> mtype;

    inline ConstSmallUpperTriMatrixViewF(const T* m, size_t s, int si, int sj) :
      mtype(m,s,si,sj) {}
    inline ConstSmallUpperTriMatrixViewF(const T* m, size_t s, int si) :
      mtype(m,s,si) {}
    inline ConstSmallUpperTriMatrixViewF(const T* m, size_t s) :
      mtype(m,s) {}
    inline ConstSmallUpperTriMatrixViewF(const type& m2) : mtype(m2) {}
    template <int Si2, int Sj2, IndexStyle I2>
    inline ConstSmallUpperTriMatrixViewF(
        const ConstUpperTriMatrixView<T,D,Si2,Sj2,C,I2>& m2) : mtype(m2) {}
    template <int Si2, int Sj2, IndexStyle I2>
    inline ConstSmallUpperTriMatrixViewF(
        const UpperTriMatrixView<T,D,Si2,Sj2,C,I2>& m2) : mtype(m2) {}
    template <int N2, int Si2, int Sj2, IndexStyle I2>
    inline ConstSmallUpperTriMatrixViewF(
        const ConstSmallUpperTriMatrixView<T,N2,D,Si2,Sj2,C,I2>& m2) :
      mtype(m2) {}
    template <int N2, int Si2, int Sj2, IndexStyle I2>
    inline ConstSmallUpperTriMatrixViewF(
        const SmallUpperTriMatrixView<T,N2,D,Si2,Sj2,C,I2>& m2) :
      mtype(m2) {}
    inline ~ConstSmallUpperTriMatrixViewF() {}

  private :
    inline void operator=(const type& m2);

  }; // ConstSmallUpperTriMatrixViewF

  template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I>
  struct Traits<SmallUpperTriMatrixView<T,N,D,Si,Sj,C,I> >
  {
    typedef T value_type;

    typedef typename Traits<T>::real_type real_type;
    typedef typename Traits<T>::complex_type complex_type;
    enum { misreal = Traits<T>::isreal };
    enum { miscomplex = Traits<T>::iscomplex };

    typedef SmallUpperTriMatrixView<T,N,D,Si,Sj,C,I> type;
    typedef const ConstSmallUpperTriMatrixView<T,N,D,Si,Sj,C,I> calc_type;
    typedef calc_type eval_type;
    typedef InvalidType inverse_type;

    enum { mcolsize = N };
    enum { mrowsize = N };
    enum { msize = N };
    enum { mshape = UpperTri };
    enum { mfort = (I == FortranStyle) };
    enum { mcalc = true };
    enum { mrowmajor = (Sj == 1) };
    enum { mcolmajor = (Si == 1) };
    enum { mstor = (mrowmajor ? RowMajor : ColMajor) };
    enum { mstepi = Si };
    enum { mstepj = Sj };
    enum { mdiagstep = IntTraits2<Si,Sj>::sum };
    enum { mconj = C };
    enum { munit = (D == UnitDiag) };

    typedef typename MCopyHelper<T,UpperTri,N,N,mrowmajor,mfort>::type copy_type;

    enum { twoSi = misreal ? Si : IntTraits<Si>::twoS };
    enum { twoSj = misreal ? Sj : IntTraits<Sj>::twoS };
    enum { notC = !C && miscomplex };

    typedef ConstVectorView<T,mstepi,C,I> const_col_type;
    typedef ConstVectorView<T,mstepi,C,I> const_col_range_type;
    typedef ConstVectorView<T,mstepj,C,I> const_row_type;
    typedef ConstVectorView<T,mstepj,C,I> const_row_range_type;
    typedef ConstVectorView<T,mdiagstep,C,I> const_diag_type;
    typedef ConstVectorView<T,mdiagstep,C,I> const_diag_range_type;

    typedef ConstUpperTriMatrixView<T,D,mstepi,mstepj,C,I> const_subtrimatrix_type;
    typedef ConstUpperTriMatrixView<T,D,UNKNOWN,UNKNOWN,C,I> const_subtrimatrix_step_type;
    typedef ConstMatrixView<T,mstepi,mstepj,C,I> const_submatrix_type;
    typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,C,I> const_submatrix_step_type;
    typedef ConstVectorView<T,UNKNOWN,C,I> const_subvector_type;

    typedef ConstSmallUpperTriMatrixView<T,N,D,mstepi,mstepj,C,I> const_view_type;
    typedef ConstSmallUpperTriMatrixView<T,N,D,mstepi,mstepj,C,CStyle> const_cview_type;
    typedef ConstSmallUpperTriMatrixView<T,N,D,mstepi,mstepj,C,FortranStyle> const_fview_type;
    typedef ConstUpperTriMatrixView<T,D,UNKNOWN,UNKNOWN,C> const_xview_type;
    typedef ConstSmallUpperTriMatrixView<T,N,D,1,mstepj,C,I> const_cmview_type;
    typedef ConstSmallUpperTriMatrixView<T,N,D,mstepi,1,C,I> const_rmview_type;
    typedef ConstSmallUpperTriMatrixView<T,N,D,mstepi,mstepj,notC,I> const_conjugate_type;
    typedef ConstSmallLowerTriMatrixView<T,N,D,mstepj,mstepi,C,I> const_transpose_type;
    typedef ConstSmallLowerTriMatrixView<T,N,D,mstepj,mstepi,notC,I> const_adjoint_type;

    typedef ConstUpperTriMatrixView<T,NonUnitDiag,mstepi,1,C,I> const_offdiag_type;
    typedef ConstSmallUpperTriMatrixView<T,N,UnitDiag,mstepi,1,C,I> const_unitdiag_type;
    typedef ConstSmallUpperTriMatrixView<real_type,N,D,twoSi,twoSj,false,I> const_realview_type;
    typedef const_realview_type const_imagview_type;
    typedef ConstSmallUpperTriMatrixView<T,N,D,mstepi,mstepj,false,I> const_nonconj_type;
    typedef MatrixView<T,mstepi,mstepj,C,I> nonconst_type;

    typedef typename AuxRef<T,C>::reference reference;

    typedef VectorView<T,mstepi,C,I> col_range_type;
    typedef VectorView<T,mstepj,C,I> row_range_type;
    typedef VectorView<T,mdiagstep,C,I> diag_type;
    typedef VectorView<T,mdiagstep,C,I> diag_range_type;

    typedef UpperTriMatrixView<T,D,mstepi,mstepj,C,I> subtrimatrix_type;
    typedef UpperTriMatrixView<T,D,UNKNOWN,UNKNOWN,C,I> subtrimatrix_step_type;
    typedef MatrixView<T,mstepi,mstepj,C,I> submatrix_type;
    typedef MatrixView<T,UNKNOWN,UNKNOWN,C,I> submatrix_step_type;
    typedef VectorView<T,UNKNOWN,C,I> subvector_type;

    typedef SmallUpperTriMatrixView<T,N,D,mstepi,mstepj,C,I> view_type;
    typedef SmallUpperTriMatrixView<T,N,D,mstepi,mstepj,C,CStyle> cview_type;
    typedef SmallUpperTriMatrixView<T,N,D,mstepi,mstepj,C,FortranStyle> fview_type;
    typedef UpperTriMatrixView<T,D,UNKNOWN,UNKNOWN,C> xview_type;
    typedef SmallUpperTriMatrixView<T,N,D,1,mstepj,C,I> cmview_type;
    typedef SmallUpperTriMatrixView<T,N,D,mstepi,1,C,I> rmview_type;
    typedef SmallUpperTriMatrixView<T,N,D,mstepi,mstepj,notC,I> conjugate_type;
    typedef SmallLowerTriMatrixView<T,N,D,mstepj,mstepi,C,I> transpose_type;
    typedef SmallLowerTriMatrixView<T,N,D,mstepj,mstepi,notC,I> adjoint_type;

    typedef UpperTriMatrixView<T,NonUnitDiag,mstepi,mstepj,C,I> offdiag_type;
    typedef SmallUpperTriMatrixView<T,N,UnitDiag,mstepi,mstepj,C,I> unitdiag_type;
    typedef SmallUpperTriMatrixView<real_type,N,D,twoSi,twoSj,false,I> realview_type;
    typedef realview_type imagview_type;
    typedef SmallUpperTriMatrixView<T,N,D,mstepi,mstepj,false,I> nonconj_type;
  };

  template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I>
  class SmallUpperTriMatrixView :
    public BaseMatrix_Tri_Mutable<SmallUpperTriMatrixView<T,N,D,Si,Sj,C,I> >
  {
  public:

    typedef SmallUpperTriMatrixView<T,N,D,Si,Sj,C,I> type;
    typedef BaseMatrix_Tri_Mutable<type> base_mut;
    typedef typename base_mut::reference reference;
    enum { mrowmajor = Traits<type>::mrowmajor };
    enum { mcolmajor = Traits<type>::mcolmajor };
    enum { mshape = Traits<type>::mshape };
    enum { munit = Traits<type>::munit };

    //
    // Constructors
    //

    inline SmallUpperTriMatrixView(T* m, size_t s, int si, int sj) :
      itsm(m), itss(s), itssi(si), itssj(sj) {}

    inline SmallUpperTriMatrixView(T* m, size_t s, int si) :
      itsm(m), itss(s), itssi(si), itssj(Sj)
    { TMVStaticAssert(Sj != UNKNOWN); }

    inline SmallUpperTriMatrixView(T* m, size_t s) :
      itsm(m), itss(s), itssi(Si), itssj(Sj)
    { TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); }

    inline SmallUpperTriMatrixView(T* m) :
      itsm(m), itss(N), itssi(Si), itssj(Sj)
    {
      TMVStaticAssert(N != UNKNOWN);
      TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
    }

    inline SmallUpperTriMatrixView(const type& m2) :
      itsm(m2.itsm), itss(m2.size()), 
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int Si2, int Sj2, IndexStyle I2>
    inline SmallUpperTriMatrixView(UpperTriMatrixView<T,D,Si2,Sj2,C,I2> m2) :
      itsm(m2.ptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int N2, int Si2, int Sj2, IndexStyle I2>
    inline SmallUpperTriMatrixView(
        SmallUpperTriMatrixView<T,N2,D,Si2,Sj2,C,I2> m2) :
      itsm(m2.ptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    inline ~SmallUpperTriMatrixView() {
#ifdef TMV_DEBUG
      itsm = 0;
#endif
    }


    //
    // Op = 
    //

    inline type& operator=(const type& m2)
    {
      base_mut::operator=(m2);
      return *this;
    }

    template <class M2>
    inline type& operator=(const BaseMatrix<M2>& m2)
    {
      base_mut::operator=(m2);
      return *this;
    }

    template <class M2>
    inline type& operator=(const BaseMatrix_Diag<M2>& m2)
    {
      base_mut::operator=(m2);
      return *this;
    }

    inline type& operator=(const T x)
    {
      base_mut::operator=(x);
      return *this;
    }


    //
    // Auxilliary Functions
    //

    inline const T* cptr() const { return itsm; }
    inline T* ptr() { return itsm; }

    inline T cref(int i, int j) const
    {
      return (
          (D==UnitDiag && i==j ) ? T(1) :
          (i>j) ? T(0) :
          DoConj<C>(itsm[i*stepi() + j*stepj()]));
    }

    inline reference ref(int i, int j)
    { return reference(itsm[i*stepi()+j*stepj()]); }

    inline size_t colsize() const { return itss; }
    inline size_t rowsize() const { return itss; }
    inline size_t size() const { return itss; }
    inline int stepi() const { return itssi; }
    inline int stepj() const { return itssj; }
    inline bool isconj() const { return C; }
    inline bool isunit() const { return D == UnitDiag; }
    inline bool isrm() const
    { return mrowmajor || (!mcolmajor &&  stepj() == 1); }
    inline bool iscm() const
    { return mcolmajor || (!mrowmajor &&  stepi() == 1); }

  private :

#ifdef TMV_DEBUG
    T* itsm;
#else
    T*const itsm;
#endif
    const StepInt<N> itss;
    const StepInt<Si> itssi;
    const StepInt<Sj> itssj;

  }; // SmallUpperTriMatrixView

  template <class T, int N, DiagType D, int Si, int Sj, bool C>
  class SmallUpperTriMatrixViewF :
    public SmallUpperTriMatrixView<T,N,D,Si,Sj,C,FortranStyle>
  {
  public:

    typedef SmallUpperTriMatrixViewF<T,N,D,Si,Sj,C> type;
    typedef SmallUpperTriMatrixView<T,N,D,Si,Sj,C,FortranStyle> mtype;

    inline SmallUpperTriMatrixViewF(T* m, size_t s, int si, int sj) :
      mtype(m,s,si,sj) {}
    inline SmallUpperTriMatrixViewF(T* m, size_t s, int si) :
      mtype(m,s,si) {}
    inline SmallUpperTriMatrixViewF(T* m, size_t s) : mtype(m,s) {}
    inline SmallUpperTriMatrixViewF(T* m) : mtype(m) {}
    inline SmallUpperTriMatrixViewF(const type& m2) : mtype(m2) {}
    template <int Si2, int Sj2, IndexStyle I2>
    inline SmallUpperTriMatrixViewF(UpperTriMatrixView<T,D,Si2,Sj2,C,I2> m2) :
      mtype(m2) {}
    template <int N2, int Si2, int Sj2, IndexStyle I2>
    inline SmallUpperTriMatrixViewF(
        SmallUpperTriMatrixView<T,N2,D,Si2,Sj2,C,I2> m2) : mtype(m2) {}
    inline ~SmallUpperTriMatrixViewF() {}

    inline type& operator=(const type& m2)
    { mtype::operator=(m2); return *this; }
    template <class M2>
    inline type& operator=(const BaseMatrix<M2>& m2)
    { mtype::operator=(m2); return *this; }
    template <class M2>
    inline type& operator=(const BaseMatrix_Diag<M2>& m2)
    { mtype::operator=(m2); return *this; }
    inline type& operator=(const T x)
    { mtype::operator=(x); return *this; }

  }; // SmallUpperTriMatrixViewF


  template <class T, int N, DiagType D, StorageType S, IndexStyle I> 
  struct Traits<SmallLowerTriMatrix<T,N,D,S,I> >
  {
    typedef T value_type;

    typedef typename Traits<T>::real_type real_type;
    typedef typename Traits<T>::complex_type complex_type;
    enum { misreal = Traits<T>::isreal };
    enum { miscomplex = Traits<T>::iscomplex };

    typedef SmallLowerTriMatrix<T,N,D,S,I> type;
    typedef const type& calc_type;
    typedef const type& eval_type;
    typedef type copy_type;

    enum { mcolsize = N };
    enum { mrowsize = N };
    enum { msize = N };
    enum { mshape = LowerTri };
    enum { mfort = (I == FortranStyle) };
    enum { mcalc = true };
    enum { mrowmajor = (S == RowMajor) };
    enum { mcolmajor = (S == ColMajor) };
    enum { mstor = S };
    enum { mstepi = (S==RowMajor ? N : 1) };
    enum { mstepj = (S==RowMajor ? 1 : N) };
    enum { mdiagstep = N+1 };
    enum { mconj = false };
    enum { munit = (D == UnitDiag) };

    enum { twoSi = misreal ? int(mstepi) : int(IntTraits<mstepi>::twoS) };
    enum { twoSj = misreal ? int(mstepj) : int(IntTraits<mstepj>::twoS) };
    enum { notC = miscomplex };

    typedef ConstVectorView<T,mstepi,false,I> const_col_range_type;
    typedef ConstVectorView<T,mstepj,false,I> const_row_range_type;
    typedef ConstSmallVectorView<T,N,mdiagstep,false,I> const_diag_type;
    typedef ConstVectorView<T,mdiagstep,false,I> const_diag_range_type;

    typedef ConstLowerTriMatrixView<T,D,mstepi,mstepj,false,I> const_subtrimatrix_type;
    typedef ConstLowerTriMatrixView<T,D,UNKNOWN,UNKNOWN,false,I> const_subtrimatrix_step_type;
    typedef ConstMatrixView<T,mstepi,mstepj,false,I> const_submatrix_type;
    typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,false,I> const_submatrix_step_type;
    typedef ConstVectorView<T,UNKNOWN,false,I> const_subvector_type;

    typedef ConstSmallLowerTriMatrixView<T,N,D,mstepi,mstepj,false,I> const_view_type;
    typedef ConstSmallLowerTriMatrixView<T,N,D,mstepi,mstepj,false,CStyle> const_cview_type;
    typedef ConstSmallLowerTriMatrixView<T,N,D,mstepi,mstepj,false,FortranStyle> const_fview_type;
    typedef ConstLowerTriMatrixView<T,D> const_xview_type;
    typedef ConstSmallLowerTriMatrixView<T,N,D,1,mstepj,false,I> const_cmview_type;
    typedef ConstSmallLowerTriMatrixView<T,N,D,mstepi,1,false,I> const_rmview_type;
    typedef ConstSmallLowerTriMatrixView<T,N,D,mstepi,mstepj,notC,I> const_conjugate_type;
    typedef ConstSmallUpperTriMatrixView<T,N,D,mstepj,mstepi,false,I> const_transpose_type;
    typedef ConstSmallUpperTriMatrixView<T,N,D,mstepj,mstepi,notC,I> const_adjoint_type;

    typedef ConstLowerTriMatrixView<T,NonUnitDiag,mstepi,mstepj,false,I> const_offdiag_type;
    typedef ConstSmallLowerTriMatrixView<T,N,UnitDiag,mstepi,mstepj,false,I> const_unitdiag_type;
    typedef ConstSmallLowerTriMatrixView<real_type,N,D,twoSi,twoSj,false,I> const_realview_type;
    typedef const_realview_type const_imagview_type;
    typedef ConstSmallLowerTriMatrixView<T,N,D,mstepi,mstepj,false,I> const_nonconj_type;

    typedef InvalidType inverse_type;
    typedef SmallLowerTriMatrixView<T,N,D,mstepi,mstepj,false,I> nonconst_type;

    typedef T& reference;

    typedef VectorView<T,mstepi,false,I> col_range_type;
    typedef VectorView<T,mstepj,false,I> row_range_type;
    typedef VectorView<T,mdiagstep,false,I> diag_type;
    typedef VectorView<T,mdiagstep,false,I> diag_range_type;

    typedef LowerTriMatrixView<T,D,mstepi,mstepj,false,I> subtrimatrix_type;
    typedef LowerTriMatrixView<T,D,UNKNOWN,UNKNOWN,false,I> subtrimatrix_step_type;
    typedef MatrixView<T,mstepi,mstepj,false,I> submatrix_type;
    typedef MatrixView<T,UNKNOWN,UNKNOWN,false,I> submatrix_step_type;
    typedef VectorView<T,UNKNOWN,false,I> subvector_type;

    typedef SmallLowerTriMatrixView<T,N,D,mstepi,mstepj,false,I> view_type;
    typedef SmallLowerTriMatrixView<T,N,D,mstepi,mstepj,false,CStyle> cview_type;
    typedef SmallLowerTriMatrixView<T,N,D,mstepi,mstepj,false,FortranStyle> fview_type;
    typedef LowerTriMatrixView<T,D> xview_type;
    typedef SmallLowerTriMatrixView<T,N,D,1,mstepj,false,I> cmview_type;
    typedef SmallLowerTriMatrixView<T,N,D,mstepi,1,false,I> rmview_type;
    typedef SmallLowerTriMatrixView<T,N,D,mstepi,mstepj,notC,I> conjugate_type;
    typedef SmallUpperTriMatrixView<T,N,D,mstepj,mstepi,false,I> transpose_type;
    typedef SmallUpperTriMatrixView<T,N,D,mstepj,mstepi,notC,I> adjoint_type;

    typedef LowerTriMatrixView<T,NonUnitDiag,mstepi,mstepj,false,I> offdiag_type;
    typedef SmallLowerTriMatrixView<T,N,UnitDiag,mstepi,mstepj,false,I> unitdiag_type;
    typedef SmallLowerTriMatrixView<real_type,N,D,twoSi,twoSj,false,I> realview_type;
    typedef realview_type imagview_type;
    typedef SmallLowerTriMatrixView<T,N,D,mstepi,mstepj,false,I> nonconj_type;
  };

#ifdef XTEST
#ifdef TMVDEBUG
#define XTEST_DEBUG
#endif
#endif

  template <class T, int N, DiagType D, StorageType S, IndexStyle I> 
  class SmallLowerTriMatrix : 
    public BaseMatrix_Tri_Mutable<SmallLowerTriMatrix<T,N,D,S,I> >
  {
  public:

    typedef SmallLowerTriMatrix<T,N,D,S,I> type;
    typedef BaseMatrix_Tri_Mutable<type> base_mut;
    enum { misreal = Traits<type>::misreal };
    enum { miscomplex = Traits<type>::miscomplex };
    enum { mshape = Traits<type>::mshape };
    enum { munit = Traits<type>::munit };
    enum { mstepi = Traits<type>::mstepi };
    enum { mstepj = Traits<type>::mstepj };

    //
    // Constructors
    //

    inline SmallLowerTriMatrix()
    {
      TMVStaticAssert(N>0);
      TMVStaticAssert(S==RowMajor || S==ColMajor); 
#ifdef TMVDEBUG
      this->SetAllTo(T(888));
#endif
    }

    explicit inline SmallLowerTriMatrix(T x)
    {
      TMVStaticAssert(N>0);
      TMVStaticAssert(S==RowMajor || S==ColMajor);
      this->SetAllTo(x);
    }

    explicit inline SmallLowerTriMatrix(const T* vv) 
    {
      TMVStaticAssert(N>0);
      TMVStaticAssert(S==RowMajor || S==ColMajor);
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      ConstSmallMatrixView<T,N,D,mstepi,mstepj> m2(vv);
      SmallMatrixView<T,N,D,mstepi,mstepj> m1(ptr());
      m2.AssignTo(m1);
    }

    explicit inline SmallLowerTriMatrix(const std::vector<T>& vv) 
    {
      TMVStaticAssert(N>0);
      TMVStaticAssert(S==RowMajor || S==ColMajor);
      TMVAssert(vv.size() == N*N);
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      std::copy(vv.begin(),vv.end(),itsm.get());
    }

    template <class M2> 
    inline SmallLowerTriMatrix(const BaseMatrix<M2>& m2) 
    { 
      TMVStaticAssert(N>0);
      TMVStaticAssert(S==RowMajor || S==ColMajor); 
      TMVStaticAssert(M2::misreal || miscomplex);
      TMVStaticAssert((ShapeTraits2<M2::mshape,mshape>::assignable));
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      m2.AssignTo(*this);
    }

    inline SmallLowerTriMatrix(const type& m2) 
    { 
      TMVStaticAssert(N>0);
      TMVStaticAssert(S==RowMajor || S==ColMajor); 
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      m2.AssignTo(*this);
    }

    template <class M2>
    inline SmallLowerTriMatrix(const BaseMatrix_Tri<M2>& m2) 
    {
      TMVStaticAssert(N>0);
      TMVStaticAssert(S==RowMajor || S==ColMajor);
      TMVStaticAssert(Shape(M2::mshape) == LowerTri);
      TMVStaticAssert(M2::munit || !munit);
      TMVStaticAssert(M2::misreal || miscomplex);
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      m2.AssignTo(*this);
      if (isunit() && !m2.isunit()) {
        if (m2.size() > 0) this->OffDiag() = m2.OffDiag();
      }
      else m2.AssignTo(*this);
    }

    template <class M2>
    inline SmallLowerTriMatrix(const BaseMatrix_Rec<M2>& m2) 
    {
      TMVStaticAssert(N>0);
      TMVStaticAssert(S==RowMajor || S==ColMajor);
      TMVStaticAssert(M2::misreal || miscomplex);
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      if (isunit()) {
        m2.UnitLowerTri().AssignTo(*this);
      } else {
        m2.LowerTri().AssignTo(*this);
      }
    }

    inline ~SmallLowerTriMatrix()
    {
#ifdef TMVDEBUG
      this->SetAllTo(T(999));
#endif
    }


    //
    // Op=
    //

    inline type& operator=(const type& m2)
    { 
      if (&m2 != this) 
        base_mut::operator=(m2);
      return *this;
    }

    template <class M2>
    inline type& operator=(const BaseMatrix<M2>& m2)
    {
      base_mut::operator=(m2);
      return *this;
    }

    template <class M2>
    inline type& operator=(const BaseMatrix_Diag<M2>& m2)
    {
      base_mut::operator=(m2);
      return *this;
    }

    inline type& operator=(T x)
    {
      base_mut::operator=(x);
      return *this;
    }


    // 
    // Auxilliary Functions
    //

    inline const T* cptr() const { return itsm; }
    inline T* ptr() { return itsm; }

    inline T cref(int i, int j) const 
    {
      return (
          (D==UnitDiag && i==j ) ? T(1) :
          (i<j) ? T(0) :
          itsm[S==RowMajor ? i*stepi() + j : i + j*stepj()]);
    }

    inline T& ref(int i, int j)
    { return itsm[S==RowMajor ? i*stepi() + j : i + j*stepj()]; }

    inline size_t size() const { return N; }
    inline int stepi() const { return mstepi; }
    inline int stepj() const { return mstepj; }
    inline DiagType dt() const { return D; }
    inline bool isunit() const { return D == UnitDiag; }
    inline bool isconj() const { return false; }
    inline bool isrm() const { return S==RowMajor; }
    inline bool iscm() const { return S==ColMajor; }
    inline StorageType stor() const { return S; }

  protected :

    StackArray<T,N*N> itsm;

  }; // SmallLowerTriMatrix

  template <class T, int N, DiagType D, StorageType S>
  class SmallLowerTriMatrixF : 
    public SmallLowerTriMatrix<T,N,D,S,FortranStyle>
  {
  public:

    typedef SmallLowerTriMatrixF<T,N,D,S> type;
    typedef SmallLowerTriMatrix<T,N,D,S,FortranStyle> mtype;

    inline SmallLowerTriMatrixF() : mtype() {}
    explicit inline SmallLowerTriMatrixF(T x) : mtype(x) {}
    explicit inline SmallLowerTriMatrixF(const T* vv) : mtype(vv) {}
    explicit inline SmallLowerTriMatrixF(const std::vector<T>& vv) :
      mtype(vv) {}
    template <class M2> 
    inline SmallLowerTriMatrixF(const BaseMatrix<M2>& m2) : mtype(m2) {}
    inline SmallLowerTriMatrixF(const type& m2) : mtype(m2) {}
    template <class M2>
    inline SmallLowerTriMatrixF(const BaseMatrix_Tri<M2>& m2) : mtype(m2) {}
    template <class M2>
    inline SmallLowerTriMatrixF(const BaseMatrix_Rec<M2>& m2) : mtype(m2) {}
    inline ~SmallLowerTriMatrixF() {}

    inline type& operator=(const type& m2)
    { mtype::operator=(m2); return *this; }
    template <class M2>
    inline type& operator=(const BaseMatrix<M2>& m2)
    { mtype::operator=(m2); return *this; }
    template <class M2>
    inline type& operator=(const BaseMatrix_Diag<M2>& m2)
    { mtype::operator=(m2); return *this; }
    inline type& operator=(T x)
    { mtype::operator=(x); return *this; }

  }; // SmallLowerTriMatrixF

  template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I>
  struct Traits<ConstSmallLowerTriMatrixView<T,N,D,Si,Sj,C,I> >
  {
    typedef T value_type;

    typedef typename Traits<T>::real_type real_type;
    typedef typename Traits<T>::complex_type complex_type;
    enum { misreal = Traits<T>::isreal };
    enum { miscomplex = Traits<T>::iscomplex };

    typedef ConstSmallLowerTriMatrixView<T,N,D,Si,Sj,C,I> type;
    typedef const type& calc_type;
    typedef const type& eval_type;
    typedef InvalidType inverse_type;

    enum { mcolsize = N };
    enum { mrowsize = N };
    enum { msize = N };
    enum { mshape = LowerTri };
    enum { mfort = (I == FortranStyle) };
    enum { mcalc = true };
    enum { mrowmajor = (Sj == 1) };
    enum { mcolmajor = (Si == 1) };
    enum { mstor = (mrowmajor ? RowMajor : ColMajor) };
    enum { mstepi = Si };
    enum { mstepj = Sj };
    enum { mdiagstep = IntTraits2<Si,Sj>::sum };
    enum { mconj = C };
    enum { munit = (D == UnitDiag) };

    // In case N == UNKNOWN
    typedef typename MCopyHelper<T,LowerTri,N,N,mrowmajor,mfort>::type copy_type;

    enum { twoSi = misreal ? Si : IntTraits<Si>::twoS };
    enum { twoSj = misreal ? Sj : IntTraits<Sj>::twoS };
    enum { notC = !C && miscomplex };

    typedef ConstVectorView<T,mstepi,C,I> const_col_range_type;
    typedef ConstVectorView<T,mstepj,C,I> const_row_range_type;
    typedef ConstVectorView<T,mdiagstep,C,I> const_diag_type;
    typedef ConstVectorView<T,mdiagstep,C,I> const_diag_range_type;

    typedef ConstLowerTriMatrixView<T,D,mstepi,mstepj,C,I> const_subtrimatrix_type;
    typedef ConstLowerTriMatrixView<T,D,UNKNOWN,UNKNOWN,C,I> const_subtrimatrix_step_type;
    typedef ConstMatrixView<T,mstepi,mstepj,C,I> const_submatrix_type;
    typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,C,I> const_submatrix_step_type;
    typedef ConstVectorView<T,UNKNOWN,C,I> const_subvector_type;

    typedef ConstSmallLowerTriMatrixView<T,N,D,mstepi,mstepj,C,I> const_view_type;
    typedef ConstSmallLowerTriMatrixView<T,N,D,mstepi,mstepj,C,CStyle> const_cview_type;
    typedef ConstSmallLowerTriMatrixView<T,N,D,mstepi,mstepj,C,FortranStyle> const_fview_type;
    typedef ConstLowerTriMatrixView<T,D,UNKNOWN,UNKNOWN,C> const_xview_type;
    typedef ConstSmallLowerTriMatrixView<T,N,D,1,mstepj,C,I> const_cmview_type;
    typedef ConstSmallLowerTriMatrixView<T,N,D,mstepi,1,C,I> const_rmview_type;
    typedef ConstSmallLowerTriMatrixView<T,N,D,mstepi,mstepj,notC,I> const_conjugate_type;
    typedef ConstSmallUpperTriMatrixView<T,N,D,mstepj,mstepi,C,I> const_transpose_type;
    typedef ConstSmallUpperTriMatrixView<T,N,D,mstepj,mstepi,notC,I> const_adjoint_type;

    typedef ConstLowerTriMatrixView<T,NonUnitDiag,mstepi,mstepj,C,I> const_offdiag_type;
    typedef ConstSmallLowerTriMatrixView<T,N,UnitDiag,mstepi,mstepj,C,I> const_unitdiag_type;
    typedef ConstSmallLowerTriMatrixView<real_type,N,D,twoSi,twoSj,false,I> const_realview_type;
    typedef const_realview_type const_imagview_type;
    typedef ConstSmallLowerTriMatrixView<T,N,D,mstepi,mstepj,false,I> const_nonconj_type;
    typedef SmallLowerTriMatrixView<T,N,D,mstepi,mstepj,C,I> nonconst_type;
  };

  template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I>
  class ConstSmallLowerTriMatrixView :
    public BaseMatrix_Tri<ConstSmallLowerTriMatrixView<T,N,D,Si,Sj,C,I> >
  {
  public:

    typedef ConstSmallLowerTriMatrixView<T,N,D,Si,Sj,C,I> type;
    enum { mrowmajor = Traits<type>::mrowmajor };
    enum { mcolmajor = Traits<type>::mcolmajor };
    enum { mshape = Traits<type>::mshape };
    enum { munit = Traits<type>::munit };

    //
    // Constructors
    //
    inline ConstSmallLowerTriMatrixView(const T* m, size_t s, int si, int sj) :
      itsm(m), itss(s), itssi(si), itssj(sj) {}

    inline ConstSmallLowerTriMatrixView(const T* m, size_t s, int si) :
      itsm(m), itss(s), itssi(si), itssj(Sj)
    { TMVStaticAssert(Sj != UNKNOWN); }

    inline ConstSmallLowerTriMatrixView(const T* m, size_t s) :
      itsm(m), itss(s), itssi(Si), itssj(Sj)
    { TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); }

    inline ConstSmallLowerTriMatrixView(const T* m) :
      itsm(m), itss(N), itssi(Si), itssj(Sj)
    {
      TMVStaticAssert(N != UNKNOWN);
      TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
    }

    inline ConstSmallLowerTriMatrixView(const type& m2) :
      itsm(m2.cptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int Si2, int Sj2, IndexStyle I2>
    inline ConstSmallLowerTriMatrixView(
        const ConstLowerTriMatrixView<T,D,Si2,Sj2,C,I2>& m2) :
      itsm(m2.cptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int Si2, int Sj2, IndexStyle I2>
    inline ConstSmallLowerTriMatrixView(
        const LowerTriMatrixView<T,D,Si2,Sj2,C,I2>& m2) :
      itsm(m2.cptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int N2, int Si2, int Sj2, IndexStyle I2>
    inline ConstSmallLowerTriMatrixView(
        const ConstSmallLowerTriMatrixView<T,N2,D,Si2,Sj2,C,I2>& m2) :
      itsm(m2.cptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int N2, int Si2, int Sj2, IndexStyle I2>
    inline ConstSmallLowerTriMatrixView(
        const SmallLowerTriMatrixView<T,N2,D,Si2,Sj2,C,I2>& m2) :
      itsm(m2.cptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    inline ~ConstSmallLowerTriMatrixView() {
#ifdef TMV_DEBUG
      itsm = 0;
#endif
    }

  private :
    inline void operator=(const type& m2);
  public :

    //
    // Auxilliary Functions
    //

    inline const T* cptr() const { return itsm; }

    inline T cref(int i, int j) const 
    {
      return (
          (D==UnitDiag && i==j ) ? T(1) :
          (i<j) ? T(0) :
          DoConj<C>(itsm[i*stepi() + j*stepj()]));
    }

    inline size_t colsize() const { return itss; }
    inline size_t rowsize() const { return itss; }
    inline size_t size() const { return itss; }
    inline int stepi() const { return itssi; }
    inline int stepj() const { return itssj; }
    inline bool isconj() const { return C; }
    inline bool isunit() const { return D == UnitDiag; }
    inline bool isrm() const
    { return mrowmajor || (!mcolmajor && stepj() == 1); }
    inline bool iscm() const
    { return mcolmajor || (!mrowmajor && stepi() == 1); }

  private :

#ifdef TMV_DEBUG
    const T* itsm;
#else
    const T*const itsm;
#endif
    const StepInt<N> itss;
    const StepInt<Si> itssi;
    const StepInt<Sj> itssj;

  }; // ConstSmallLowerTriMatrixView

  template <class T, int N, DiagType D, int Si, int Sj, bool C>
  class ConstSmallLowerTriMatrixViewF :
    public ConstSmallLowerTriMatrixView<T,N,D,Si,Sj,C,FortranStyle>
  {
  public:

    typedef ConstSmallLowerTriMatrixViewF<T,N,D,Si,Sj,C> type;
    typedef ConstSmallLowerTriMatrixView<T,N,D,Si,Sj,C,FortranStyle> mtype;

    inline ConstSmallLowerTriMatrixViewF(const T* m, size_t s, int si, int sj) :
      mtype(m,s,si,sj) {}
    inline ConstSmallLowerTriMatrixViewF(const T* m, size_t s, int si) :
      mtype(m,s,si) {}
    inline ConstSmallLowerTriMatrixViewF(const T* m, size_t s) :
      mtype(m,s) {}
    inline ConstSmallLowerTriMatrixViewF(const type& m2) : mtype(m2) {}
    template <int Si2, int Sj2, IndexStyle I2>
    inline ConstSmallLowerTriMatrixViewF(
        const ConstLowerTriMatrixView<T,D,Si2,Sj2,C,I2>& m2) : mtype(m2) {}
    template <int Si2, int Sj2, IndexStyle I2>
    inline ConstSmallLowerTriMatrixViewF(
        const LowerTriMatrixView<T,D,Si2,Sj2,C,I2>& m2) : mtype(m2) {}
    template <int N2, int Si2, int Sj2, IndexStyle I2>
    inline ConstSmallLowerTriMatrixViewF(
        const ConstSmallLowerTriMatrixView<T,N2,D,Si2,Sj2,C,I2>& m2) :
      mtype(m2) {}
    template <int N2, int Si2, int Sj2, IndexStyle I2>
    inline ConstSmallLowerTriMatrixViewF(
        const SmallLowerTriMatrixView<T,N2,D,Si2,Sj2,C,I2>& m2) :
      mtype(m2) {}
    inline ~ConstSmallLowerTriMatrixViewF() {}

  private :
    inline void operator=(const type& m2);

  }; // ConstSmallLowerTriMatrixViewF

  template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I>
  struct Traits<SmallLowerTriMatrixView<T,N,D,Si,Sj,C,I> >
  {
    typedef T value_type;

    typedef typename Traits<T>::real_type real_type;
    typedef typename Traits<T>::complex_type complex_type;
    enum { misreal = Traits<T>::isreal };
    enum { miscomplex = Traits<T>::iscomplex };

    typedef SmallLowerTriMatrixView<T,N,D,Si,Sj,C,I> type;
    typedef const ConstSmallLowerTriMatrixView<T,N,D,Si,Sj,C,I> calc_type;
    typedef calc_type eval_type;
    typedef InvalidType inverse_type;

    enum { mcolsize = N };
    enum { mrowsize = N };
    enum { msize = N };
    enum { mshape = LowerTri };
    enum { mfort = (I == FortranStyle) };
    enum { mcalc = true };
    enum { mrowmajor = (Sj == 1) };
    enum { mcolmajor = (Si == 1) };
    enum { mstor = (mrowmajor ? RowMajor : ColMajor) };
    enum { mstepi = Si };
    enum { mstepj = Sj };
    enum { mdiagstep = IntTraits2<Si,Sj>::sum };
    enum { mconj = C };
    enum { munit = (D == UnitDiag) };

    typedef typename MCopyHelper<T,LowerTri,N,N,mrowmajor,mfort>::type copy_type;

    enum { twoSi = misreal ? Si : IntTraits<Si>::twoS };
    enum { twoSj = misreal ? Sj : IntTraits<Sj>::twoS };
    enum { notC = !C && miscomplex };

    typedef ConstVectorView<T,mstepi,C,I> const_col_type;
    typedef ConstVectorView<T,mstepi,C,I> const_col_range_type;
    typedef ConstVectorView<T,mstepj,C,I> const_row_type;
    typedef ConstVectorView<T,mstepj,C,I> const_row_range_type;
    typedef ConstVectorView<T,mdiagstep,C,I> const_diag_type;
    typedef ConstVectorView<T,mdiagstep,C,I> const_diag_range_type;

    typedef ConstLowerTriMatrixView<T,D,mstepi,mstepj,C,I> const_subtrimatrix_type;
    typedef ConstLowerTriMatrixView<T,D,UNKNOWN,UNKNOWN,C,I> const_subtrimatrix_step_type;
    typedef ConstMatrixView<T,mstepi,mstepj,C,I> const_submatrix_type;
    typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,C,I> const_submatrix_step_type;
    typedef ConstVectorView<T,UNKNOWN,C,I> const_subvector_type;

    typedef ConstSmallLowerTriMatrixView<T,N,D,mstepi,mstepj,C,I> const_view_type;
    typedef ConstSmallLowerTriMatrixView<T,N,D,mstepi,mstepj,C,CStyle> const_cview_type;
    typedef ConstSmallLowerTriMatrixView<T,N,D,mstepi,mstepj,C,FortranStyle> const_fview_type;
    typedef ConstLowerTriMatrixView<T,D,UNKNOWN,UNKNOWN,C> const_xview_type;
    typedef ConstSmallLowerTriMatrixView<T,N,D,1,mstepj,C,I> const_cmview_type;
    typedef ConstSmallLowerTriMatrixView<T,N,D,mstepi,1,C,I> const_rmview_type;
    typedef ConstSmallLowerTriMatrixView<T,N,D,mstepi,mstepj,notC,I> const_conjugate_type;
    typedef ConstSmallUpperTriMatrixView<T,N,D,mstepj,mstepi,C,I> const_transpose_type;
    typedef ConstSmallUpperTriMatrixView<T,N,D,mstepj,mstepi,notC,I> const_adjoint_type;

    typedef ConstLowerTriMatrixView<T,NonUnitDiag,mstepi,1,C,I> const_offdiag_type;
    typedef ConstSmallLowerTriMatrixView<T,N,UnitDiag,mstepi,1,C,I> const_unitdiag_type;
    typedef ConstSmallLowerTriMatrixView<real_type,N,D,twoSi,twoSj,false,I> const_realview_type;
    typedef const_realview_type const_imagview_type;
    typedef ConstSmallLowerTriMatrixView<T,N,D,mstepi,mstepj,false,I> const_nonconj_type;
    typedef MatrixView<T,mstepi,mstepj,C,I> nonconst_type;

    typedef typename AuxRef<T,C>::reference reference;

    typedef VectorView<T,mstepi,C,I> col_range_type;
    typedef VectorView<T,mstepj,C,I> row_range_type;
    typedef VectorView<T,mdiagstep,C,I> diag_type;
    typedef VectorView<T,mdiagstep,C,I> diag_range_type;

    typedef LowerTriMatrixView<T,D,mstepi,mstepj,C,I> subtrimatrix_type;
    typedef LowerTriMatrixView<T,D,UNKNOWN,UNKNOWN,C,I> subtrimatrix_step_type;
    typedef MatrixView<T,mstepi,mstepj,C,I> submatrix_type;
    typedef MatrixView<T,UNKNOWN,UNKNOWN,C,I> submatrix_step_type;
    typedef VectorView<T,UNKNOWN,C,I> subvector_type;

    typedef SmallLowerTriMatrixView<T,N,D,mstepi,mstepj,C,I> view_type;
    typedef SmallLowerTriMatrixView<T,N,D,mstepi,mstepj,C,CStyle> cview_type;
    typedef SmallLowerTriMatrixView<T,N,D,mstepi,mstepj,C,FortranStyle> fview_type;
    typedef LowerTriMatrixView<T,D,UNKNOWN,UNKNOWN,C> xview_type;
    typedef SmallLowerTriMatrixView<T,N,D,1,mstepj,C,I> cmview_type;
    typedef SmallLowerTriMatrixView<T,N,D,mstepi,1,C,I> rmview_type;
    typedef SmallLowerTriMatrixView<T,N,D,mstepi,mstepj,notC,I> conjugate_type;
    typedef SmallUpperTriMatrixView<T,N,D,mstepj,mstepi,C,I> transpose_type;
    typedef SmallUpperTriMatrixView<T,N,D,mstepj,mstepi,notC,I> adjoint_type;

    typedef LowerTriMatrixView<T,NonUnitDiag,mstepi,mstepj,C,I> offdiag_type;
    typedef SmallLowerTriMatrixView<T,N,UnitDiag,mstepi,mstepj,C,I> unitdiag_type;
    typedef SmallLowerTriMatrixView<real_type,N,D,twoSi,twoSj,false,I> realview_type;
    typedef realview_type imagview_type;
    typedef SmallLowerTriMatrixView<T,N,D,mstepi,mstepj,false,I> nonconj_type;
  };

  template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I>
  class SmallLowerTriMatrixView :
    public BaseMatrix_Tri_Mutable<SmallLowerTriMatrixView<T,N,D,Si,Sj,C,I> >
  {
  public:

    typedef SmallLowerTriMatrixView<T,N,D,Si,Sj,C,I> type;
    typedef BaseMatrix_Tri_Mutable<type> base_mut;
    typedef typename base_mut::reference reference;
    enum { mrowmajor = Traits<type>::mrowmajor };
    enum { mcolmajor = Traits<type>::mcolmajor };
    enum { mshape = Traits<type>::mshape };
    enum { munit = Traits<type>::munit };

    //
    // Constructors
    //

    inline SmallLowerTriMatrixView(T* m, size_t s, int si, int sj) :
      itsm(m), itss(s), itssi(si), itssj(sj) {}

    inline SmallLowerTriMatrixView(T* m, size_t s, int si) :
      itsm(m), itss(s), itssi(si), itssj(Sj)
    { TMVStaticAssert(Sj != UNKNOWN); }

    inline SmallLowerTriMatrixView(T* m, size_t s) :
      itsm(m), itss(s), itssi(Si), itssj(Sj)
    { TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); }

    inline SmallLowerTriMatrixView(T* m) :
      itsm(m), itss(N), itssi(Si), itssj(Sj)
    {
      TMVStaticAssert(N != UNKNOWN);
      TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
    }

    inline SmallLowerTriMatrixView(const type& m2) :
      itsm(m2.itsm), itss(m2.size()), 
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int Si2, int Sj2, IndexStyle I2>
    inline SmallLowerTriMatrixView(LowerTriMatrixView<T,D,Si2,Sj2,C,I2> m2) :
      itsm(m2.ptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int N2, int Si2, int Sj2, IndexStyle I2>
    inline SmallLowerTriMatrixView(
        SmallLowerTriMatrixView<T,N2,D,Si2,Sj2,C,I2> m2) :
      itsm(m2.ptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    inline ~SmallLowerTriMatrixView() {
#ifdef TMV_DEBUG
      itsm = 0;
#endif
    }


    //
    // Op = 
    //

    inline type& operator=(const type& m2)
    {
      base_mut::operator=(m2);
      return *this;
    }

    template <class M2>
    inline type& operator=(const BaseMatrix<M2>& m2)
    {
      base_mut::operator=(m2);
      return *this;
    }

    template <class M2>
    inline type& operator=(const BaseMatrix_Diag<M2>& m2)
    {
      base_mut::operator=(m2);
      return *this;
    }

    inline type& operator=(const T x)
    {
      base_mut::operator=(x);
      return *this;
    }


    //
    // Auxilliary Functions
    //

    inline const T* cptr() const { return itsm; }
    inline T* ptr() { return itsm; }

    inline T cref(int i, int j) const
    {
      return (
          (D==UnitDiag && i==j ) ? T(1) :
          (i<j) ? T(0) :
          DoConj<C>(itsm[i*stepi() + j*stepj()]));
    }

    inline reference ref(int i, int j)
    { return reference(itsm[i*stepi()+j*stepj()]); }

    inline size_t colsize() const { return itss; }
    inline size_t rowsize() const { return itss; }
    inline size_t size() const { return itss; }
    inline int stepi() const { return itssi; }
    inline int stepj() const { return itssj; }
    inline bool isconj() const { return C; }
    inline bool isunit() const { return D == UnitDiag; }
    inline bool isrm() const
    { return mrowmajor || (!mcolmajor &&  stepj() == 1); }
    inline bool iscm() const
    { return mcolmajor || (!mrowmajor &&  stepi() == 1); }

  private :

#ifdef TMV_DEBUG
    T* itsm;
#else
    T*const itsm;
#endif
    const StepInt<N> itss;
    const StepInt<Si> itssi;
    const StepInt<Sj> itssj;

  }; // SmallLowerTriMatrixView

  template <class T, int N, DiagType D, int Si, int Sj, bool C>
  class SmallLowerTriMatrixViewF :
    public SmallLowerTriMatrixView<T,N,D,Si,Sj,C,FortranStyle>
  {
  public:

    typedef SmallLowerTriMatrixViewF<T,N,D,Si,Sj,C> type;
    typedef SmallLowerTriMatrixView<T,N,D,Si,Sj,C,FortranStyle> mtype;

    inline SmallLowerTriMatrixViewF(T* m, size_t s, int si, int sj) :
      mtype(m,s,si,sj) {}
    inline SmallLowerTriMatrixViewF(T* m, size_t s, int si) :
      mtype(m,s,si) {}
    inline SmallLowerTriMatrixViewF(T* m, size_t s) : mtype(m,s) {}
    inline SmallLowerTriMatrixViewF(T* m) : mtype(m) {}
    inline SmallLowerTriMatrixViewF(const type& m2) : mtype(m2) {}
    template <int Si2, int Sj2, IndexStyle I2>
    inline SmallLowerTriMatrixViewF(LowerTriMatrixView<T,D,Si2,Sj2,C,I2> m2) :
      mtype(m2) {}
    template <int N2, int Si2, int Sj2, IndexStyle I2>
    inline SmallLowerTriMatrixViewF(
        SmallLowerTriMatrixView<T,N2,D,Si2,Sj2,C,I2> m2) : mtype(m2) {}
    inline ~SmallLowerTriMatrixViewF() {}

    inline type& operator=(const type& m2)
    { mtype::operator=(m2); return *this; }
    template <class M2>
    inline type& operator=(const BaseMatrix<M2>& m2)
    { mtype::operator=(m2); return *this; }
    template <class M2>
    inline type& operator=(const BaseMatrix_Diag<M2>& m2)
    { mtype::operator=(m2); return *this; }
    inline type& operator=(const T x)
    { mtype::operator=(x); return *this; }

  }; // SmallLowerTriMatrixViewF


  //
  // Swap
  //

  template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I, class M>
  inline void Swap(
      BaseMatrix_Tri_Mutable<M>& m1,
      SmallUpperTriMatrixView<T,N,D,Si,Sj,C,I> m2)
  { DoSwap(m1,m2); }
  template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I, class M>
  inline void Swap(
      SmallUpperTriMatrixView<T,N,D,Si,Sj,C,I> m1,
      BaseMatrix_Tri_Mutable<M>& m2)
  { DoSwap(m1,m2); }
  template <class T, int N, DiagType D, int Si1, int Sj1, bool C1, IndexStyle I1, int Si2, int Sj2, bool C2, IndexStyle I2>
  inline void Swap(
      SmallUpperTriMatrixView<T,N,D,Si1,Sj1,C1,I1> m1,
      SmallUpperTriMatrixView<T,N,D,Si2,Sj2,C2,I2> m2)
  { DoSwap(m1,m2); }
  template <class T, int N, DiagType D, int Si1, int Sj1, bool C1, IndexStyle I1, int Si2, int Sj2, bool C2, IndexStyle I2>
  inline void Swap(
      SmallUpperTriMatrixView<T,N,D,Si1,Sj1,C1,I1> m1,
      UpperTriMatrixView<T,D,Si2,Sj2,C2,I2> m2)
  { DoSwap(m1,m2); }
  template <class T, int N, DiagType D, int Si1, int Sj1, bool C1, IndexStyle I1, int Si2, int Sj2, bool C2, IndexStyle I2>
  inline void Swap(
      UpperTriMatrixView<T,D,Si1,Sj1,C1,I1> m1,
      SmallUpperTriMatrixView<T,N,D,Si2,Sj2,C2,I2> m2)
  { DoSwap(m1,m2); }

  template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I, class M>
  inline void Swap(
      BaseMatrix_Tri_Mutable<M>& m1,
      SmallLowerTriMatrixView<T,N,D,Si,Sj,C,I> m2)
  { Swap(m1.Transpose(),m2.Transpose()); }
  template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I, class M>
  inline void Swap(
      SmallLowerTriMatrixView<T,N,D,Si,Sj,C,I> m1,
      BaseMatrix_Tri_Mutable<M>& m2)
  { Swap(m1.Transpose(),m2.Transpose()); }
  template <class T, int N, DiagType D, int Si1, int Sj1, bool C1, IndexStyle I1, int Si2, int Sj2, bool C2, IndexStyle I2>
  inline void Swap(
      SmallLowerTriMatrixView<T,N,D,Si1,Sj1,C1,I1> m1,
      SmallLowerTriMatrixView<T,N,D,Si2,Sj2,C2,I2> m2)
  { Swap(m1.Transpose(),m2.Transpose()); }
  template <class T, int N, DiagType D, int Si1, int Sj1, bool C1, IndexStyle I1, int Si2, int Sj2, bool C2, IndexStyle I2>
  inline void Swap(
      SmallLowerTriMatrixView<T,N,D,Si1,Sj1,C1,I1> m1,
      LowerTriMatrixView<T,D,Si2,Sj2,C2,I2> m2)
  { Swap(m1.Transpose(),m2.Transpose()); }
  template <class T, int N, DiagType D, int Si1, int Sj1, bool C1, IndexStyle I1, int Si2, int Sj2, bool C2, IndexStyle I2>
  inline void Swap(
      LowerTriMatrixView<T,D,Si1,Sj1,C1,I1> m1,
      SmallLowerTriMatrixView<T,N,D,Si2,Sj2,C2,I2> m2)
  { Swap(m1.Transpose(),m2.Transpose()); }


  //
  // TypeText 
  //

  template <class T, int N, DiagType D, StorageType S, IndexStyle I>
  inline std::string TypeText(const SmallUpperTriMatrix<T,N,D,S,I>& )
  {
    std::ostringstream s;
    s << "SmallUpperTriMatrix<"<<TypeText(T())<<","<<N<<","<<Text(D);
    s << ","<<Text(S)<<","<<Text(I)<<">";
    return s.str();
  }

  template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I>
  inline std::string TypeText(const ConstSmallUpperTriMatrixView<T,N,D,Si,Sj,C,I>& m)
  {
    std::ostringstream s;
    s << "ConstSmallUpperTriMatrixView<"<<TypeText(T());
    s << ","<<IntTraits<N>::text();
    if (N == UNKNOWN) s << "("<<m.stepi()<<")";
    s <<","<<Text(D);
    s << ","<<IntTraits<Si>::text();
    if (Si == UNKNOWN) s << "("<<m.stepi()<<")";
    s << ","<<IntTraits<Sj>::text();
    if (Sj == UNKNOWN) s << "("<<m.stepj()<<")";
    s << ","<<C<<","<<Text(I)<<">";
    return s.str();
  }

  template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I>
  inline std::string TypeText(const SmallUpperTriMatrixView<T,N,D,Si,Sj,C,I>& m)
  {
    std::ostringstream s;
    s << "SmallUpperTriMatrixView<"<<TypeText(T());
    s << ","<<IntTraits<N>::text();
    if (N == UNKNOWN) s << "("<<m.stepi()<<")";
    s <<","<<Text(D);
    s << ","<<IntTraits<Si>::text();
    if (Si == UNKNOWN) s << "("<<m.stepi()<<")";
    s << ","<<IntTraits<Sj>::text();
    if (Sj == UNKNOWN) s << "("<<m.stepj()<<")";
    s << ","<<C<<","<<Text(I)<<">";
    return s.str();
  }

  template <class T, int N, DiagType D, StorageType S, IndexStyle I>
  inline std::string TypeText(const SmallLowerTriMatrix<T,N,D,S,I>& )
  {
    std::ostringstream s;
    s << "SmallLowerTriMatrix<"<<TypeText(T())<<","<<N<<","<<Text(D);
    s << ","<<Text(S)<<","<<Text(I)<<">";
    return s.str();
  }

  template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I>
  inline std::string TypeText(const ConstSmallLowerTriMatrixView<T,N,D,Si,Sj,C,I>& m)
  {
    std::ostringstream s;
    s << "ConstSmallLowerTriMatrixView<"<<TypeText(T());
    s << ","<<IntTraits<N>::text();
    if (N == UNKNOWN) s << "("<<m.stepi()<<")";
    s <<","<<Text(D);
    s << ","<<IntTraits<Si>::text();
    if (Si == UNKNOWN) s << "("<<m.stepi()<<")";
    s << ","<<IntTraits<Sj>::text();
    if (Sj == UNKNOWN) s << "("<<m.stepj()<<")";
    s << ","<<C<<","<<Text(I)<<">";
    return s.str();
  }

  template <class T, int N, DiagType D, int Si, int Sj, bool C, IndexStyle I>
  inline std::string TypeText(const SmallLowerTriMatrixView<T,N,D,Si,Sj,C,I>& m)
  {
    std::ostringstream s;
    s << "SmallLowerTriMatrixView<"<<TypeText(T());
    s << ","<<IntTraits<N>::text();
    if (N == UNKNOWN) s << "("<<m.stepi()<<")";
    s <<","<<Text(D);
    s << ","<<IntTraits<Si>::text();
    if (Si == UNKNOWN) s << "("<<m.stepi()<<")";
    s << ","<<IntTraits<Sj>::text();
    if (Sj == UNKNOWN) s << "("<<m.stepj()<<")";
    s << ","<<C<<","<<Text(I)<<">";
    return s.str();
  }


} // namespace tmv

#endif
