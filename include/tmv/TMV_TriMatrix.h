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
// This file defines the TMV TriMatrix class.
//
// Constructors:
//
//    There are two TriMatrix classes: UpperTriMatrix<T> and LowerTriMatrix<T>
//    For these notes, I will just write TriMatrix, but for all uses,
//    you need to write "Upper" or "Lower" before the "Tri".
//
//    In addition to the type template parameter (T), TriMatrixes have three
//    additional template parameters:
//        DiagType dt = UnitDiag || NonUnitDiag 
//        StorageType stor = RowMajor || ColMajor
//        IndexType I = CStyle || FortranStyle
//
//        They each have default values, so you can omit all three,
//        stor and I, or just stor.  
//        The default values are: {NonUnitDiag, RowMajor, CStyle}
//
//        If dt is UnitDiag, then the diagonal elements are not
//        actually stored or referenced.  The are all taken to be = 1.
//
//        The storage follows the same meaning as for regular Matrices.
//
//    TriMatrix<T,dt,stor,I>(size_t n)
//        Makes a Triangular Matrix with column size = row size = n
//        with _uninitialized_ values.
//
//    TriMatrix<T,dt,stor,I>(size_t n, T x)
//        Makes a Triangular Matrix with column size = row size = n
//        with all values = x
//
//    TriMatrix<T,dt,stor,I>(size_t n, T* vv)
//    TriMatrix<T,dt,stor,I>(size_t n, const std::vector<T>& vv)
//        Makes a Triangular Matrix with column size = row size = n
//        which copies the values from vv.
//
//    TriMatrix<T,dt,stor,I>(const Matrix<T>& m)
//    TriMatrix<T,dt,stor,I>(const TriMatrix<T>& m)
//        Makes a TriMatrix which copies the corresponding elements of m.
//
//
// Special Creators:
//
//    ConstUpperTriMatrixView UpperTriMatrixViewOf(
//            const T* m, size_t size, StorageType stor)
//    ConstUpperTriMatrixView UnitUpperTriMatrixViewOf(
//            const T* m, size_t size, StorageType stor)
//    UpperTriMatrixView UnitUpperTriMatrixViewOf(
//            T* m, size_t size, StorageType stor)
//    UpperTriMatrixView UnitUnitUpperTriMatrixViewOf(
//            T* m, size_t size, StorageType stor)
//        Returns a TriMatrixView of the elements in m, using the 
//        actual elements m for the storage.  The Unit versions return
//        views with dt = UnitDiag, the non-Unit versions return views
//        with dt = NonUnitDiag.
//        There are also corresponding LowerTriMatrix versions of these.
//        
//
// Access Functions
//
//    size_t colsize() const
//    size_t rowsize() const
//    size_t size() const
//        Return the dimensions of the TriMatrix
//
//    value_type operator()(int i, int j) const
//    value_type cref(int i, int j) const
//        Return the (i,j) element of the TriMatrix
//        The first one respects the index-style of the underlying matrix.
//        The second, cref, always uses CStyle indexing and does not 
//        do any checking of the valididty of i,j.
//
//    reference operator()(int i, int j)
//    reference ref(int i, int j)
//        Return a reference to the (i,j) element of the matrix
//        The first one respects the index-style of the underlying matrix.
//        The second, ref, always uses CStyle indexing and does not 
//        do any checking of the valididty of i,j.
//   
//    row_range_type row(int i, int j1, int j2)
//    const_row_range_type row(int i, int j1, int j2) const
//        Return a portion of the ith row 
//        This range must be a valid range for the requested row.
//
//    col_range_type col(int j, int i1, int i2)
//    const_col_range_type col(int j, int i1, int i2) const
//        Return a portion of the jth column
//        This range must be a valid range for the requested column.
//
//    diag_type diag()
//    const_diag_type diag() const
//        Return the main diagonal
//        The TriMatrix must be NonUnitDiag.
//
//    diag_range_type diag(int i)
//    diag_range_type diag(int i, int j1, int j2)
//    const_diag_range_type diag(int i) const
//    const_diag_range_type diag(int i, int j1, int j2) const
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
// Functions of Matrices:
//
//    Most of these are the same as for a regular matrix, so I 
//    only give the full description for the new functionality.
//
//    value_type m.Det() const    or Det(m)
//    real_type m.LogDet(value_type* sign=0) const   or LogDet(m,sign)
//    value_type m.Trace() const    or Trace(m)
//    real_type m.Norm() const    or Norm(m)
//    real_type m.NormF() const    or NormF(m)
//    real_type m.NormSq() const    or NormSq()
//    real_type m.NormSq(real_type scale) const
//    real_type m.Norm1() const    or Norm1(m)
//    real_type m.Norm2() const    or Norm2(m)
//    real_type m.NormInf() const    or NormInf(m)
//    value_type SumElements() const    or SumElements(m) 
//    real_type SumAbsElements() const    or SumAbsElements(m) 
//    value_type MaxElement() const    or MaxElement(m) 
//    value_type MinElement() const    or MinElement(m) 
//    real_type MaxAbsElement() const    or MaxAbsElement(m) 
//    real_type MinAbsElement() const    or MinAbsElement(m) 
//
//    void m.Inverse(minv) const
//        This function allows minv to be either a regular Matrix
//        or a TriMatrix (of the same Upper or Lower as m).
//    void m.InverseATA(invata) const
//    inverse_type m.Inverse() const    or Inverse(m)
//
//
// Modifying Functions
//
//    type& Zero()
//    type& SetAllTo(value_type x)
//    type& AddToAll(value_type x)
//    type& Clip(real_type thresh)
//    type& ApplyToAll(const F& f)
//    type& ConjugateSelf()
//    type& SetToIdentity(value_type x = 1)
//    Swap(TriMatrix& m1, TriMatrix& m2)
//
//    type& InvertSelf()
//        Change the TriMatrix into its own inverse.  This can be done
//        efficiently in place without requiring extra storage, so 
//        this function provides that functionality.
//
//
// Views of a TriMatrix:
//
//    (As usual, all of these have a const_version as well.)
//
//    submatrix_type SubMatrix(int i1, int i2, int j1, int j2,
//            int istep=1, int jstep=1)
//        This member function will return a submatrix using rows i1 to i2
//        and columns j1 to j2 which refers
//        to the same physical elements as the original.
//        The submatrix must be completely contained within the TriMatrix.
//
//    subvector_type SubVector(int i, int j, int istep, int jstep, int size)
//        Returns a SubVector which starts at position (i,j) in the 
//        matrix, moves in the directions (istep,jstep) and has a length
//        of size.
//
//    subtrimatrix_type SubTriMatrix(int i1, int i2, int istep)
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
//    offdiag_type OffDiag()
//        Returns the (NonUnitDiag) TriMatrix of all the off-diagonal
//        elements of a TriMatrix.
//
//    unitdiag_type ViewAsUnitDiag()
//        Re-view a NonUnitDiag TriMatrix with a view that takes the 
//        diagonal elements to all be equal to 1.
//
//    realview_type Real()
//    imagview_type Imag()
//        For a complex TriMatrix, returns the real or imaginary part
//        as a real TriMatrix.
//
//    view_type View()
//    conjugate_type Conjugate()
//    transpose_type Transpose()
//    adjoint_type Adjoint()
//        Note that the Transpose or Adjoint of an UpperTriMatrix returns 
//        a view which is a LowerTriMatrix, and vice versa.
//
//    nonconj_type NonConj()
//        Returns a mutable view of the (const) TriMatrix
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

#include "tmv/TMV_BaseMatrix_Tri.h"
#include "tmv/TMV_Matrix.h"
#include <vector>

namespace tmv {

  template <class T, DiagType D, StorageType S, IndexStyle I> 
  struct Traits<UpperTriMatrix<T,D,S,I> >
  {
    typedef T value_type;

    typedef typename Traits<T>::real_type real_type;
    typedef typename Traits<T>::complex_type complex_type;
    enum { misreal = Traits<T>::isreal };
    enum { miscomplex = Traits<T>::iscomplex };

    typedef UpperTriMatrix<T,D,S,I> type;
    typedef const type& calc_type;
    typedef const type& eval_type;
    typedef type copy_type;

    enum { mcolsize = UNKNOWN };
    enum { mrowsize = UNKNOWN };
    enum { msize = UNKNOWN };
    enum { mfort = (I == FortranStyle) };
    enum { mcalc = true };
    enum { mrowmajor = (S == RowMajor) };
    enum { mcolmajor = (S == ColMajor) };
    enum { mstor = S };
    enum { mstepi = (S==ColMajor ? 1 : UNKNOWN) };
    enum { mstepj = (S==RowMajor ? 1 : UNKNOWN) };
    enum { mdiagstep = UNKNOWN };
    enum { mconj = false };
    enum { munit = (D == UnitDiag) };
    enum { mshape = munit ? UnitUpperTri : UpperTri };

    enum { twoSi = misreal ? int(mstepi) : int(IntTraits<mstepi>::twoS) };
    enum { twoSj = misreal ? int(mstepj) : int(IntTraits<mstepj>::twoS) };
    enum { notC = miscomplex };

    typedef ConstVectorView<T,mstepi,false,I> const_col_range_type;
    typedef ConstVectorView<T,mstepj,false,I> const_row_range_type;
    typedef ConstVectorView<T,mdiagstep,false,I> const_diag_type;
    typedef ConstVectorView<T,mdiagstep,false,I> const_diag_range_type;

    typedef ConstUpperTriMatrixView<T,D,mstepi,mstepj,false,I> const_subtrimatrix_type;
    typedef ConstUpperTriMatrixView<T,D,UNKNOWN,UNKNOWN,false,I> const_subtrimatrix_step_type;
    typedef ConstMatrixView<T,mstepi,mstepj,false,I> const_submatrix_type;
    typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,false,I> const_submatrix_step_type;
    typedef ConstVectorView<T,UNKNOWN,false,I> const_subvector_type;

    typedef ConstUpperTriMatrixView<T,D,mstepi,mstepj,false,I> const_view_type;
    typedef ConstUpperTriMatrixView<T,D,mstepi,mstepj,false,CStyle> const_cview_type;
    typedef ConstUpperTriMatrixView<T,D,mstepi,mstepj,false,FortranStyle> const_fview_type;
    typedef ConstUpperTriMatrixView<T,D> const_xview_type;
    typedef ConstUpperTriMatrixView<T,D,1,mstepj,false,I> const_cmview_type;
    typedef ConstUpperTriMatrixView<T,D,mstepi,1,false,I> const_rmview_type;
    typedef ConstUpperTriMatrixView<T,D,mstepi,mstepj,notC,I> const_conjugate_type;
    typedef ConstLowerTriMatrixView<T,D,mstepj,mstepi,false,I> const_transpose_type;
    typedef ConstLowerTriMatrixView<T,D,mstepj,mstepi,notC,I> const_adjoint_type;

    typedef ConstUpperTriMatrixView<T,NonUnitDiag,mstepi,mstepj,false,I> const_offdiag_type;
    typedef ConstUpperTriMatrixView<T,UnitDiag,mstepi,mstepj,false,I> const_unitdiag_type;
    typedef ConstUpperTriMatrixView<real_type,D,twoSi,twoSj,false,I> const_realview_type;
    typedef const_realview_type const_imagview_type;
    typedef ConstUpperTriMatrixView<T,D,mstepi,mstepj,false,I> const_nonconj_type;

    typedef InvalidType inverse_type;
    typedef UpperTriMatrixView<T,D,mstepi,mstepj,false,I> nonconst_type;

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

    typedef UpperTriMatrixView<T,D,mstepi,mstepj,false,I> view_type;
    typedef UpperTriMatrixView<T,D,mstepi,mstepj,false,CStyle> cview_type;
    typedef UpperTriMatrixView<T,D,mstepi,mstepj,false,FortranStyle> fview_type;
    typedef UpperTriMatrixView<T,D> xview_type;
    typedef UpperTriMatrixView<T,D,1,mstepj,false,I> cmview_type;
    typedef UpperTriMatrixView<T,D,mstepi,1,false,I> rmview_type;
    typedef UpperTriMatrixView<T,D,mstepi,mstepj,notC,I> conjugate_type;
    typedef LowerTriMatrixView<T,D,mstepj,mstepi,false,I> transpose_type;
    typedef LowerTriMatrixView<T,D,mstepj,mstepi,notC,I> adjoint_type;

    typedef UpperTriMatrixView<T,NonUnitDiag,mstepi,mstepj,false,I> offdiag_type;
    typedef UpperTriMatrixView<T,UnitDiag,mstepi,mstepj,false,I> unitdiag_type;
    typedef UpperTriMatrixView<real_type,D,twoSi,twoSj,false,I> realview_type;
    typedef realview_type imagview_type;
    typedef UpperTriMatrixView<T,D,mstepi,mstepj,false,I> nonconj_type;
  };

#ifdef XTEST
#ifdef TMVDEBUG
#define XTEST_DEBUG
#endif
#endif

  template <class T, DiagType D, StorageType S, IndexStyle I> 
  class UpperTriMatrix : 
    public BaseMatrix_Tri_Mutable<UpperTriMatrix<T,D,S,I> >
  {
  public:

    typedef UpperTriMatrix<T,D,S,I> type;
    typedef BaseMatrix_Tri_Mutable<type> base_mut;
    enum { misreal = Traits<type>::misreal };
    enum { miscomplex = Traits<type>::miscomplex };
    enum { mshape = Traits<type>::mshape };
    enum { munit = Traits<type>::munit };

    //
    // Constructors
    //

    explicit inline UpperTriMatrix(size_t _size) :
      itss(_size), itsm(new T[itss*itss])
    {
      TMVStaticAssert(S==RowMajor || S==ColMajor); 
#ifdef TMVDEBUG
      this->SetAllTo(T(888));
#endif
    }

    inline UpperTriMatrix(size_t _size, T x) : 
      itss(_size), itsm(new T[itss*itss])
    {
      TMVStaticAssert(S==RowMajor || S==ColMajor);
      this->SetAllTo(x);
    }

    inline UpperTriMatrix(size_t _size, const T* vv) :
      itss(_size), itsm(new T[itss*itss])
    {
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      TMVStaticAssert(S==RowMajor || S==ColMajor);
      std::copy(vv,vv+itss*itss,itsm.get());
    }

    inline UpperTriMatrix(size_t _size, const std::vector<T>& vv) :
      itss(_size), itsm(new T[itss*itss])
    {
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      TMVStaticAssert(S==RowMajor || S==ColMajor);
      TMVAssert(vv.size() == itss*itss);
      std::copy(vv.begin(),vv.end(),itsm.get());
    }

    template <class M2> 
    inline UpperTriMatrix(const BaseMatrix<M2>& m2) :
      itss(m2.rowsize()), itsm(new T[itss*itss])
    { 
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      TMVStaticAssert(S==RowMajor || S==ColMajor); 
      TMVStaticAssert(M2::misreal || miscomplex);
      TMVStaticAssert((ShapeTraits2<M2::mshape,mshape>::assignable));
      m2.AssignTo(*this);
    }

    inline UpperTriMatrix(const type& m2) :
      itss(m2.size()), itsm(new T[itss*itss])
    { 
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      TMVStaticAssert(S==RowMajor || S==ColMajor); 
      m2.AssignTo(*this);
    }

    template <class M2>
    inline UpperTriMatrix(const BaseMatrix_Tri<M2>& m2) :
      itss(m2.size()), itsm(new T[itss*itss])
    {
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      TMVStaticAssert(S==RowMajor || S==ColMajor);
      TMVStaticAssert((ShapeTraits<M2::mshape>::upper));
      TMVStaticAssert(M2::misreal || miscomplex);
      Maybe<munit>::offdiag_copy(m2,*this);
    }

    template <class M2>
    inline UpperTriMatrix(const BaseMatrix_Rec<M2>& m2) :
      itss(m2.rowsize()), itsm(new T[itss*itss])
    {
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      TMVStaticAssert(S==RowMajor || S==ColMajor);
      TMVStaticAssert(M2::misreal || miscomplex);
      Maybe<D==UnitDiag>::unit_uppertri(m2).AssignTo(*this);
    }

    inline ~UpperTriMatrix()
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

    inline const T* cptr() const { return itsm.get(); }
    inline T* ptr() { return itsm.get(); }

    inline T cref(int i, int j) const 
    {
      return (
          (D==UnitDiag && i==j ) ? T(1) :
          (i>j) ? T(0) :
          itsm[S==RowMajor ? i*stepi() + j : i + j*stepj()]);
    }

    inline T& ref(int i, int j)
    { return itsm[S==RowMajor ? i*stepi() + j : i + j*stepj()]; }

    inline size_t size() const { return itss; }
    inline int stepi() const { return S==RowMajor ? itss : 1; }
    inline int stepj() const { return S==RowMajor ? 1 : itss; }
    inline DiagType dt() const { return D; }
    inline bool isunit() const { return D == UnitDiag; }
    inline bool isconj() const { return false; }
    inline bool isrm() const { return S==RowMajor; }
    inline bool iscm() const { return S==ColMajor; }
    inline StorageType stor() const { return S; }

  protected :

    const size_t itss;
    auto_array<T> itsm;

  }; // UpperTriMatrix

  template <class T, DiagType D, StorageType S>
  class UpperTriMatrixF : 
    public UpperTriMatrix<T,D,S,FortranStyle>
  {
  public:

    typedef UpperTriMatrixF<T,D,S> type;
    typedef UpperTriMatrix<T,D,S,FortranStyle> mtype;

    explicit inline UpperTriMatrixF(size_t s) : mtype(s) {}
    inline UpperTriMatrixF(size_t s, T x) : mtype(s,x) {}
    inline UpperTriMatrixF(size_t s, const T* vv) : mtype(s,vv) {}
    inline UpperTriMatrixF(size_t s, const std::vector<T>& vv) : mtype(s,vv) {}
    template <class M2> 
    inline UpperTriMatrixF(const BaseMatrix<M2>& m2) : mtype(m2) {}
    inline UpperTriMatrixF(const type& m2) : mtype(m2) {}
    template <class M2>
    inline UpperTriMatrixF(const BaseMatrix_Tri<M2>& m2) : mtype(m2) {}
    template <class M2>
    inline UpperTriMatrixF(const BaseMatrix_Rec<M2>& m2) : mtype(m2) {}
    inline ~UpperTriMatrixF() {}

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

  }; // UpperTriMatrixF

  template <class T, DiagType D, int Si, int Sj, bool C, IndexStyle I>
  struct Traits<ConstUpperTriMatrixView<T,D,Si,Sj,C,I> >
  {
    typedef T value_type;

    typedef typename Traits<T>::real_type real_type;
    typedef typename Traits<T>::complex_type complex_type;
    enum { misreal = Traits<T>::isreal };
    enum { miscomplex = Traits<T>::iscomplex };

    typedef ConstUpperTriMatrixView<T,D,Si,Sj,C,I> type;
    typedef const type& calc_type;
    typedef const type& eval_type;
    typedef UpperTriMatrix<T,D,Sj==1?RowMajor:ColMajor,I> copy_type;
    typedef InvalidType inverse_type;

    enum { mcolsize = UNKNOWN };
    enum { mrowsize = UNKNOWN };
    enum { msize = UNKNOWN };
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
    enum { mshape = munit ? UnitUpperTri : UpperTri };

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

    typedef ConstUpperTriMatrixView<T,D,mstepi,mstepj,C,I> const_view_type;
    typedef ConstUpperTriMatrixView<T,D,mstepi,mstepj,C,CStyle> const_cview_type;
    typedef ConstUpperTriMatrixView<T,D,mstepi,mstepj,C,FortranStyle> const_fview_type;
    typedef ConstUpperTriMatrixView<T,D,UNKNOWN,UNKNOWN,C> const_xview_type;
    typedef ConstUpperTriMatrixView<T,D,1,mstepj,C,I> const_cmview_type;
    typedef ConstUpperTriMatrixView<T,D,mstepi,1,C,I> const_rmview_type;
    typedef ConstUpperTriMatrixView<T,D,mstepi,mstepj,notC,I> const_conjugate_type;
    typedef ConstLowerTriMatrixView<T,D,mstepj,mstepi,C,I> const_transpose_type;
    typedef ConstLowerTriMatrixView<T,D,mstepj,mstepi,notC,I> const_adjoint_type;

    typedef ConstUpperTriMatrixView<T,NonUnitDiag,mstepi,mstepj,C,I> const_offdiag_type;
    typedef ConstUpperTriMatrixView<T,UnitDiag,mstepi,mstepj,C,I> const_unitdiag_type;
    typedef ConstUpperTriMatrixView<real_type,D,twoSi,twoSj,false,I> const_realview_type;
    typedef const_realview_type const_imagview_type;
    typedef ConstUpperTriMatrixView<T,D,mstepi,mstepj,false,I> const_nonconj_type;
    typedef UpperTriMatrixView<T,D,mstepi,mstepj,C,I> nonconst_type;
  };

  template <class T, DiagType D, int Si, int Sj, bool C, IndexStyle I>
  class ConstUpperTriMatrixView :
    public BaseMatrix_Tri<ConstUpperTriMatrixView<T,D,Si,Sj,C,I> >
  {
  public:

    typedef ConstUpperTriMatrixView<T,D,Si,Sj,C,I> type;
    enum { mrowmajor = Traits<type>::mrowmajor };
    enum { mcolmajor = Traits<type>::mcolmajor };
    enum { mshape = Traits<type>::mshape };
    enum { munit = Traits<type>::munit };

    //
    // Constructors
    //
    inline ConstUpperTriMatrixView(const T* m, size_t s, int si, int sj) :
      itsm(m), itss(s), itssi(si), itssj(sj) {}

    inline ConstUpperTriMatrixView(const T* m, size_t s, int si) :
      itsm(m), itss(s), itssi(si), itssj(Sj)
    { TMVStaticAssert(Sj != UNKNOWN); }

    inline ConstUpperTriMatrixView(const T* m, size_t s) :
      itsm(m), itss(s), itssi(Si), itssj(Sj)
    { TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); }

    inline ConstUpperTriMatrixView(const type& m2) :
      itsm(m2.cptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int Si2, int Sj2, IndexStyle I2>
    inline ConstUpperTriMatrixView(
        const ConstUpperTriMatrixView<T,D,Si2,Sj2,C,I2>& m2) :
      itsm(m2.cptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int Si2, int Sj2, IndexStyle I2>
    inline ConstUpperTriMatrixView(
        const UpperTriMatrixView<T,D,Si2,Sj2,C,I2>& m2) :
      itsm(m2.cptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int N2, int Si2, int Sj2, IndexStyle I2>
    inline ConstUpperTriMatrixView(
        const ConstSmallUpperTriMatrixView<T,N2,D,Si2,Sj2,C,I2>& m2) :
      itsm(m2.cptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int N2, int Si2, int Sj2, IndexStyle I2>
    inline ConstUpperTriMatrixView(
        const SmallUpperTriMatrixView<T,N2,D,Si2,Sj2,C,I2>& m2) :
      itsm(m2.cptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    inline ~ConstUpperTriMatrixView() {
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
    const size_t itss;
    const StepInt<Si> itssi;
    const StepInt<Sj> itssj;

  }; // ConstUpperTriMatrixView

  template <class T, DiagType D, int Si, int Sj, bool C>
  class ConstUpperTriMatrixViewF :
    public ConstUpperTriMatrixView<T,D,Si,Sj,C,FortranStyle>
  {
  public:

    typedef ConstUpperTriMatrixViewF<T,D,Si,Sj,C> type;
    typedef ConstUpperTriMatrixView<T,D,Si,Sj,C,FortranStyle> mtype;

    inline ConstUpperTriMatrixViewF(const T* m, size_t s, int si, int sj) :
      mtype(m,s,si,sj) {}
    inline ConstUpperTriMatrixViewF(const T* m, size_t s, int si) :
      mtype(m,s,si) {}
    inline ConstUpperTriMatrixViewF(const T* m, size_t s) :
      mtype(m,s) {}
    inline ConstUpperTriMatrixViewF(const type& m2) : mtype(m2) {}
    template <int Si2, int Sj2, IndexStyle I2>
    inline ConstUpperTriMatrixViewF(
        const ConstUpperTriMatrixView<T,D,Si2,Sj2,C,I2>& m2) : mtype(m2) {}
    template <int Si2, int Sj2, IndexStyle I2>
    inline ConstUpperTriMatrixViewF(
        const UpperTriMatrixView<T,D,Si2,Sj2,C,I2>& m2) : mtype(m2) {}
    template <int N2, int Si2, int Sj2, IndexStyle I2>
    inline ConstUpperTriMatrixViewF(
        const ConstSmallUpperTriMatrixView<T,N2,D,Si2,Sj2,C,I2>& m2) :
      mtype(m2) {}
    template <int N2, int Si2, int Sj2, IndexStyle I2>
    inline ConstUpperTriMatrixViewF(
        const SmallUpperTriMatrixView<T,N2,D,Si2,Sj2,C,I2>& m2) :
      mtype(m2) {}
    inline ~ConstUpperTriMatrixViewF() {}

  private :
    inline void operator=(const type& m2);

  }; // ConstUpperTriMatrixViewF

  template <class T, DiagType D, int Si, int Sj, bool C, IndexStyle I>
  struct Traits<UpperTriMatrixView<T,D,Si,Sj,C,I> >
  {
    typedef T value_type;

    typedef typename Traits<T>::real_type real_type;
    typedef typename Traits<T>::complex_type complex_type;
    enum { misreal = Traits<T>::isreal };
    enum { miscomplex = Traits<T>::iscomplex };

    typedef UpperTriMatrixView<T,D,Si,Sj,C,I> type;
    typedef const ConstUpperTriMatrixView<T,D,Si,Sj,C,I> calc_type;
    typedef calc_type eval_type;
    typedef UpperTriMatrix<T,D,Sj==1?RowMajor:ColMajor,I> copy_type;
    typedef InvalidType inverse_type;

    enum { mcolsize = UNKNOWN };
    enum { mrowsize = UNKNOWN };
    enum { msize = UNKNOWN };
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
    enum { mshape = munit ? UnitUpperTri : UpperTri };

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

    typedef ConstUpperTriMatrixView<T,D,mstepi,mstepj,C,I> const_view_type;
    typedef ConstUpperTriMatrixView<T,D,mstepi,mstepj,C,CStyle> const_cview_type;
    typedef ConstUpperTriMatrixView<T,D,mstepi,mstepj,C,FortranStyle> const_fview_type;
    typedef ConstUpperTriMatrixView<T,D,UNKNOWN,UNKNOWN,C> const_xview_type;
    typedef ConstUpperTriMatrixView<T,D,1,mstepj,C,I> const_cmview_type;
    typedef ConstUpperTriMatrixView<T,D,mstepi,1,C,I> const_rmview_type;
    typedef ConstUpperTriMatrixView<T,D,mstepi,mstepj,notC,I> const_conjugate_type;
    typedef ConstLowerTriMatrixView<T,D,mstepj,mstepi,C,I> const_transpose_type;
    typedef ConstLowerTriMatrixView<T,D,mstepj,mstepi,notC,I> const_adjoint_type;

    typedef ConstUpperTriMatrixView<T,NonUnitDiag,mstepi,1,C,I> const_offdiag_type;
    typedef ConstUpperTriMatrixView<T,UnitDiag,mstepi,1,C,I> const_unitdiag_type;
    typedef ConstUpperTriMatrixView<real_type,D,twoSi,twoSj,false,I> const_realview_type;
    typedef const_realview_type const_imagview_type;
    typedef ConstUpperTriMatrixView<T,D,mstepi,mstepj,false,I> const_nonconj_type;
    typedef UpperTriMatrixView<T,D,mstepi,mstepj,C,I> nonconst_type;

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

    typedef UpperTriMatrixView<T,D,mstepi,mstepj,C,I> view_type;
    typedef UpperTriMatrixView<T,D,mstepi,mstepj,C,CStyle> cview_type;
    typedef UpperTriMatrixView<T,D,mstepi,mstepj,C,FortranStyle> fview_type;
    typedef UpperTriMatrixView<T,D,UNKNOWN,UNKNOWN,C> xview_type;
    typedef UpperTriMatrixView<T,D,1,mstepj,C,I> cmview_type;
    typedef UpperTriMatrixView<T,D,mstepi,1,C,I> rmview_type;
    typedef UpperTriMatrixView<T,D,mstepi,mstepj,notC,I> conjugate_type;
    typedef LowerTriMatrixView<T,D,mstepj,mstepi,C,I> transpose_type;
    typedef LowerTriMatrixView<T,D,mstepj,mstepi,notC,I> adjoint_type;

    typedef UpperTriMatrixView<T,NonUnitDiag,mstepi,mstepj,C,I> offdiag_type;
    typedef UpperTriMatrixView<T,UnitDiag,mstepi,mstepj,C,I> unitdiag_type;
    typedef UpperTriMatrixView<real_type,D,twoSi,twoSj,false,I> realview_type;
    typedef realview_type imagview_type;
    typedef UpperTriMatrixView<T,D,mstepi,mstepj,false,I> nonconj_type;
  };

  template <class T, DiagType D, int Si, int Sj, bool C, IndexStyle I>
  class UpperTriMatrixView :
    public BaseMatrix_Tri_Mutable<UpperTriMatrixView<T,D,Si,Sj,C,I> >
  {
  public:

    typedef UpperTriMatrixView<T,D,Si,Sj,C,I> type;
    typedef BaseMatrix_Tri_Mutable<type> base_mut;
    typedef typename base_mut::reference reference;
    enum { mrowmajor = Traits<type>::mrowmajor };
    enum { mcolmajor = Traits<type>::mcolmajor };
    enum { mshape = Traits<type>::mshape };
    enum { munit = Traits<type>::munit };

    //
    // Constructors
    //

    inline UpperTriMatrixView(T* m, size_t s, int si, int sj) :
      itsm(m), itss(s), itssi(si), itssj(sj) {}

    inline UpperTriMatrixView(T* m, size_t s, int si) :
      itsm(m), itss(s), itssi(si), itssj(Sj)
    { TMVStaticAssert(Sj != UNKNOWN); }

    inline UpperTriMatrixView(T* m, size_t s) :
      itsm(m), itss(s), itssi(Si), itssj(Sj)
    { TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); }

    inline UpperTriMatrixView(const type& m2) :
      itsm(m2.itsm), itss(m2.size()), 
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int Si2, int Sj2, IndexStyle I2>
    inline UpperTriMatrixView(UpperTriMatrixView<T,D,Si2,Sj2,C,I2> m2) :
      itsm(m2.ptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int N2, int Si2, int Sj2, IndexStyle I2>
    inline UpperTriMatrixView(SmallUpperTriMatrixView<T,N2,D,Si2,Sj2,C,I2> m2) :
      itsm(m2.ptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    inline ~UpperTriMatrixView() {
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
    const size_t itss;
    const StepInt<Si> itssi;
    const StepInt<Sj> itssj;

  }; // UpperTriMatrixView

  template <class T, DiagType D, int Si, int Sj, bool C>
  class UpperTriMatrixViewF :
    public UpperTriMatrixView<T,D,Si,Sj,C,FortranStyle>
  {
  public:

    typedef UpperTriMatrixViewF<T,D,Si,Sj,C> type;
    typedef UpperTriMatrixView<T,D,Si,Sj,C,FortranStyle> mtype;

    inline UpperTriMatrixViewF(T* m, size_t s, int si, int sj) :
      mtype(m,s,si,sj) {}
    inline UpperTriMatrixViewF(T* m, size_t s, int si) :
      mtype(m,s,si) {}
    inline UpperTriMatrixViewF(T* m, size_t s) : mtype(m,s) {}
    inline UpperTriMatrixViewF(const type& m2) : mtype(m2) {}
    template <int Si2, int Sj2, IndexStyle I2>
    inline UpperTriMatrixViewF(UpperTriMatrixView<T,D,Si2,Sj2,C,I2> m2) :
      mtype(m2) {}
    template <int N2, int Si2, int Sj2, IndexStyle I2>
    inline UpperTriMatrixViewF(
        SmallUpperTriMatrixView<T,N2,D,Si2,Sj2,C,I2> m2) : mtype(m2) {}
    inline ~UpperTriMatrixViewF() {}

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

  }; // UpperTriMatrixViewF


  template <class T, DiagType D, StorageType S, IndexStyle I> 
  struct Traits<LowerTriMatrix<T,D,S,I> >
  {
    typedef T value_type;

    typedef typename Traits<T>::real_type real_type;
    typedef typename Traits<T>::complex_type complex_type;
    enum { misreal = Traits<T>::isreal };
    enum { miscomplex = Traits<T>::iscomplex };

    typedef LowerTriMatrix<T,D,S,I> type;
    typedef const type& calc_type;
    typedef const type& eval_type;
    typedef type copy_type;

    enum { mcolsize = UNKNOWN };
    enum { mrowsize = UNKNOWN };
    enum { msize = UNKNOWN };
    enum { mfort = (I == FortranStyle) };
    enum { mcalc = true };
    enum { mrowmajor = (S == RowMajor) };
    enum { mcolmajor = (S == ColMajor) };
    enum { mstor = S };
    enum { mstepi = (S==ColMajor ? 1 : UNKNOWN) };
    enum { mstepj = (S==RowMajor ? 1 : UNKNOWN) };
    enum { mdiagstep = UNKNOWN };
    enum { mconj = false };
    enum { munit = (D == UnitDiag) };
    enum { mshape = munit ? UnitLowerTri : LowerTri };

    enum { twoSi = misreal ? int(mstepi) : int(IntTraits<mstepi>::twoS) };
    enum { twoSj = misreal ? int(mstepj) : int(IntTraits<mstepj>::twoS) };
    enum { notC = miscomplex };

    typedef ConstVectorView<T,mstepi,false,I> const_col_range_type;
    typedef ConstVectorView<T,mstepj,false,I> const_row_range_type;
    typedef ConstVectorView<T,mdiagstep,false,I> const_diag_type;
    typedef ConstVectorView<T,mdiagstep,false,I> const_diag_range_type;

    typedef ConstLowerTriMatrixView<T,D,mstepi,mstepj,false,I> const_subtrimatrix_type;
    typedef ConstLowerTriMatrixView<T,D,UNKNOWN,UNKNOWN,false,I> const_subtrimatrix_step_type;
    typedef ConstMatrixView<T,mstepi,mstepj,false,I> const_submatrix_type;
    typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,false,I> const_submatrix_step_type;
    typedef ConstVectorView<T,UNKNOWN,false,I> const_subvector_type;

    typedef ConstLowerTriMatrixView<T,D,mstepi,mstepj,false,I> const_view_type;
    typedef ConstLowerTriMatrixView<T,D,mstepi,mstepj,false,CStyle> const_cview_type;
    typedef ConstLowerTriMatrixView<T,D,mstepi,mstepj,false,FortranStyle> const_fview_type;
    typedef ConstLowerTriMatrixView<T,D> const_xview_type;
    typedef ConstLowerTriMatrixView<T,D,1,mstepj,false,I> const_cmview_type;
    typedef ConstLowerTriMatrixView<T,D,mstepi,1,false,I> const_rmview_type;
    typedef ConstLowerTriMatrixView<T,D,mstepi,mstepj,notC,I> const_conjugate_type;
    typedef ConstUpperTriMatrixView<T,D,mstepj,mstepi,false,I> const_transpose_type;
    typedef ConstUpperTriMatrixView<T,D,mstepj,mstepi,notC,I> const_adjoint_type;

    typedef ConstLowerTriMatrixView<T,NonUnitDiag,mstepi,mstepj,false,I> const_offdiag_type;
    typedef ConstLowerTriMatrixView<T,UnitDiag,mstepi,mstepj,false,I> const_unitdiag_type;
    typedef ConstLowerTriMatrixView<real_type,D,twoSi,twoSj,false,I> const_realview_type;
    typedef const_realview_type const_imagview_type;
    typedef ConstLowerTriMatrixView<T,D,mstepi,mstepj,false,I> const_nonconj_type;

    typedef InvalidType inverse_type;
    typedef LowerTriMatrixView<T,D,mstepi,mstepj,false,I> nonconst_type;

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

    typedef LowerTriMatrixView<T,D,mstepi,mstepj,false,I> view_type;
    typedef LowerTriMatrixView<T,D,mstepi,mstepj,false,CStyle> cview_type;
    typedef LowerTriMatrixView<T,D,mstepi,mstepj,false,FortranStyle> fview_type;
    typedef LowerTriMatrixView<T,D> xview_type;
    typedef LowerTriMatrixView<T,D,1,mstepj,false,I> cmview_type;
    typedef LowerTriMatrixView<T,D,mstepi,1,false,I> rmview_type;
    typedef LowerTriMatrixView<T,D,mstepi,mstepj,notC,I> conjugate_type;
    typedef UpperTriMatrixView<T,D,mstepj,mstepi,false,I> transpose_type;
    typedef UpperTriMatrixView<T,D,mstepj,mstepi,notC,I> adjoint_type;

    typedef LowerTriMatrixView<T,NonUnitDiag,mstepi,mstepj,false,I> offdiag_type;
    typedef LowerTriMatrixView<T,UnitDiag,mstepi,mstepj,false,I> unitdiag_type;
    typedef LowerTriMatrixView<real_type,D,twoSi,twoSj,false,I> realview_type;
    typedef realview_type imagview_type;
    typedef LowerTriMatrixView<T,D,mstepi,mstepj,false,I> nonconj_type;
  };

#ifdef XTEST
#ifdef TMVDEBUG
#define XTEST_DEBUG
#endif
#endif

  template <class T, DiagType D, StorageType S, IndexStyle I> 
  class LowerTriMatrix : 
    public BaseMatrix_Tri_Mutable<LowerTriMatrix<T,D,S,I> >
  {
  public:

    typedef LowerTriMatrix<T,D,S,I> type;
    typedef BaseMatrix_Tri_Mutable<type> base_mut;
    enum { misreal = Traits<type>::misreal };
    enum { miscomplex = Traits<type>::miscomplex };
    enum { mshape = Traits<type>::mshape };
    enum { munit = Traits<type>::munit };

    //
    // Constructors
    //

    explicit inline LowerTriMatrix(size_t _size) :
      itss(_size), itsm(new T[itss*itss])
    {
      TMVStaticAssert(S==RowMajor || S==ColMajor); 
#ifdef TMVDEBUG
      this->SetAllTo(T(888));
#endif
    }

    inline LowerTriMatrix(size_t _size, T x) : 
      itss(_size), itsm(new T[itss*itss])
    {
      TMVStaticAssert(S==RowMajor || S==ColMajor);
      this->SetAllTo(x);
    }

    inline LowerTriMatrix(size_t _size, const T* vv) :
      itss(_size), itsm(new T[itss*itss])
    {
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      TMVStaticAssert(S==RowMajor || S==ColMajor);
      std::copy(vv,vv+itss*itss,itsm.get());
    }

    inline LowerTriMatrix(size_t _size, const std::vector<T>& vv) :
      itss(_size), itsm(new T[itss*itss])
    {
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      TMVStaticAssert(S==RowMajor || S==ColMajor);
      TMVAssert(vv.size() == itss*itss);
      std::copy(vv.begin(),vv.end(),itsm.get());
    }

    template <class M2> 
    inline LowerTriMatrix(const BaseMatrix<M2>& m2) :
      itss(m2.rowsize()), itsm(new T[itss*itss])
    { 
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      TMVStaticAssert(S==RowMajor || S==ColMajor); 
      TMVStaticAssert(M2::misreal || miscomplex);
      TMVStaticAssert((ShapeTraits2<M2::mshape,mshape>::assignable));
      m2.AssignTo(*this);
    }

    inline LowerTriMatrix(const type& m2) :
      itss(m2.size()), itsm(new T[itss*itss])
    { 
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      TMVStaticAssert(S==RowMajor || S==ColMajor); 
      m2.AssignTo(*this);
    }

    template <class M2>
    inline LowerTriMatrix(const BaseMatrix_Tri<M2>& m2) :
      itss(m2.size()), itsm(new T[itss*itss])
    {
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      TMVStaticAssert(S==RowMajor || S==ColMajor);
      TMVStaticAssert((ShapeTraits<M2::mshape>::lower));
      TMVStaticAssert(M2::misreal || miscomplex);
      Maybe<munit>::offdiag_copy(m2,*this);
    }

    template <class M2>
    inline LowerTriMatrix(const BaseMatrix_Rec<M2>& m2) :
      itss(m2.colsize()), itsm(new T[itss*itss])
    {
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      TMVStaticAssert(S==RowMajor || S==ColMajor);
      TMVStaticAssert(M2::misreal || miscomplex);
      Maybe<D==UnitDiag>::unit_lowertri(m2).AssignTo(*this);
    }

    inline ~LowerTriMatrix()
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

    inline const T* cptr() const { return itsm.get(); }
    inline T* ptr() { return itsm.get(); }

    inline T cref(int i, int j) const 
    {
      return (
          (D==UnitDiag && i==j ) ? T(1) :
          (i<j) ? T(0) :
          itsm[S==RowMajor ? i*stepi() + j : i + j*stepj()]);
    }

    inline T& ref(int i, int j)
    { return itsm[S==RowMajor ? i*stepi() + j : i + j*stepj()]; }

    inline size_t size() const { return itss; }
    inline int stepi() const { return S==RowMajor ? itss : 1; }
    inline int stepj() const { return S==RowMajor ? 1 : itss; }
    inline DiagType dt() const { return D; }
    inline bool isunit() const { return D == UnitDiag; }
    inline bool isconj() const { return false; }
    inline bool isrm() const { return S==RowMajor; }
    inline bool iscm() const { return S==ColMajor; }
    inline StorageType stor() const { return S; }

  protected :

    const size_t itss;
    auto_array<T> itsm;

  }; // LowerTriMatrix

  template <class T, DiagType D, StorageType S>
  class LowerTriMatrixF : 
    public LowerTriMatrix<T,D,S,FortranStyle>
  {
  public:

    typedef LowerTriMatrixF<T,D,S> type;
    typedef LowerTriMatrix<T,D,S,FortranStyle> mtype;

    explicit inline LowerTriMatrixF(size_t s) : mtype(s) {}
    inline LowerTriMatrixF(size_t s, T x) : mtype(s,x) {}
    inline LowerTriMatrixF(size_t s, const T* vv) : mtype(s,vv) {}
    inline LowerTriMatrixF(size_t s, const std::vector<T>& vv) : mtype(s,vv) {}
    template <class M2> 
    inline LowerTriMatrixF(const BaseMatrix<M2>& m2) : mtype(m2) {}
    inline LowerTriMatrixF(const type& m2) : mtype(m2) {}
    template <class M2>
    inline LowerTriMatrixF(const BaseMatrix_Tri<M2>& m2) : mtype(m2) {}
    template <class M2>
    inline LowerTriMatrixF(const BaseMatrix_Rec<M2>& m2) : mtype(m2) {}
    inline ~LowerTriMatrixF() {}

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

  }; // LowerTriMatrixF

  template <class T, DiagType D, int Si, int Sj, bool C, IndexStyle I>
  struct Traits<ConstLowerTriMatrixView<T,D,Si,Sj,C,I> >
  {
    typedef T value_type;

    typedef typename Traits<T>::real_type real_type;
    typedef typename Traits<T>::complex_type complex_type;
    enum { misreal = Traits<T>::isreal };
    enum { miscomplex = Traits<T>::iscomplex };

    typedef ConstLowerTriMatrixView<T,D,Si,Sj,C,I> type;
    typedef const type& calc_type;
    typedef const type& eval_type;
    typedef LowerTriMatrix<T,D,Sj==1?RowMajor:ColMajor,I> copy_type;
    typedef InvalidType inverse_type;

    enum { mcolsize = UNKNOWN };
    enum { mrowsize = UNKNOWN };
    enum { msize = UNKNOWN };
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
    enum { mshape = munit ? UnitLowerTri : LowerTri };

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

    typedef ConstLowerTriMatrixView<T,D,mstepi,mstepj,C,I> const_view_type;
    typedef ConstLowerTriMatrixView<T,D,mstepi,mstepj,C,CStyle> const_cview_type;
    typedef ConstLowerTriMatrixView<T,D,mstepi,mstepj,C,FortranStyle> const_fview_type;
    typedef ConstLowerTriMatrixView<T,D,UNKNOWN,UNKNOWN,C> const_xview_type;
    typedef ConstLowerTriMatrixView<T,D,1,mstepj,C,I> const_cmview_type;
    typedef ConstLowerTriMatrixView<T,D,mstepi,1,C,I> const_rmview_type;
    typedef ConstLowerTriMatrixView<T,D,mstepi,mstepj,notC,I> const_conjugate_type;
    typedef ConstUpperTriMatrixView<T,D,mstepj,mstepi,C,I> const_transpose_type;
    typedef ConstUpperTriMatrixView<T,D,mstepj,mstepi,notC,I> const_adjoint_type;

    typedef ConstLowerTriMatrixView<T,NonUnitDiag,mstepi,mstepj,C,I> const_offdiag_type;
    typedef ConstLowerTriMatrixView<T,UnitDiag,mstepi,mstepj,C,I> const_unitdiag_type;
    typedef ConstLowerTriMatrixView<real_type,D,twoSi,twoSj,false,I> const_realview_type;
    typedef const_realview_type const_imagview_type;
    typedef ConstLowerTriMatrixView<T,D,mstepi,mstepj,false,I> const_nonconj_type;
    typedef LowerTriMatrixView<T,D,mstepi,mstepj,C,I> nonconst_type;
  };

  template <class T, DiagType D, int Si, int Sj, bool C, IndexStyle I>
  class ConstLowerTriMatrixView :
    public BaseMatrix_Tri<ConstLowerTriMatrixView<T,D,Si,Sj,C,I> >
  {
  public:

    typedef ConstLowerTriMatrixView<T,D,Si,Sj,C,I> type;
    enum { mrowmajor = Traits<type>::mrowmajor };
    enum { mcolmajor = Traits<type>::mcolmajor };
    enum { mshape = Traits<type>::mshape };
    enum { munit = Traits<type>::munit };

    //
    // Constructors
    //
    inline ConstLowerTriMatrixView(const T* m, size_t s, int si, int sj) :
      itsm(m), itss(s), itssi(si), itssj(sj) {}

    inline ConstLowerTriMatrixView(const T* m, size_t s, int si) :
      itsm(m), itss(s), itssi(si), itssj(Sj)
    { TMVStaticAssert(Sj != UNKNOWN); }

    inline ConstLowerTriMatrixView(const T* m, size_t s) :
      itsm(m), itss(s), itssi(Si), itssj(Sj)
    { TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); }

    inline ConstLowerTriMatrixView(const type& m2) :
      itsm(m2.cptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int Si2, int Sj2, IndexStyle I2>
    inline ConstLowerTriMatrixView(
        const ConstLowerTriMatrixView<T,D,Si2,Sj2,C,I2>& m2) :
      itsm(m2.cptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int Si2, int Sj2, IndexStyle I2>
    inline ConstLowerTriMatrixView(
        const LowerTriMatrixView<T,D,Si2,Sj2,C,I2>& m2) :
      itsm(m2.cptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int N2, int Si2, int Sj2, IndexStyle I2>
    inline ConstLowerTriMatrixView(
        const ConstSmallLowerTriMatrixView<T,N2,D,Si2,Sj2,C,I2>& m2) :
      itsm(m2.cptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int N2, int Si2, int Sj2, IndexStyle I2>
    inline ConstLowerTriMatrixView(
        const SmallLowerTriMatrixView<T,N2,D,Si2,Sj2,C,I2>& m2) :
      itsm(m2.cptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    inline ~ConstLowerTriMatrixView() {
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
    const size_t itss;
    const StepInt<Si> itssi;
    const StepInt<Sj> itssj;

  }; // ConstLowerTriMatrixView

  template <class T, DiagType D, int Si, int Sj, bool C>
  class ConstLowerTriMatrixViewF :
    public ConstLowerTriMatrixView<T,D,Si,Sj,C,FortranStyle>
  {
  public:

    typedef ConstLowerTriMatrixViewF<T,D,Si,Sj,C> type;
    typedef ConstLowerTriMatrixView<T,D,Si,Sj,C,FortranStyle> mtype;

    inline ConstLowerTriMatrixViewF(const T* m, size_t s, int si, int sj) :
      mtype(m,s,si,sj) {}
    inline ConstLowerTriMatrixViewF(const T* m, size_t s, int si) :
      mtype(m,s,si) {}
    inline ConstLowerTriMatrixViewF(const T* m, size_t s) :
      mtype(m,s) {}
    inline ConstLowerTriMatrixViewF(const type& m2) : mtype(m2) {}
    template <int Si2, int Sj2, IndexStyle I2>
    inline ConstLowerTriMatrixViewF(
        const ConstLowerTriMatrixView<T,D,Si2,Sj2,C,I2>& m2) : mtype(m2) {}
    template <int Si2, int Sj2, IndexStyle I2>
    inline ConstLowerTriMatrixViewF(
        const LowerTriMatrixView<T,D,Si2,Sj2,C,I2>& m2) : mtype(m2) {}
    template <int N2, int Si2, int Sj2, IndexStyle I2>
    inline ConstLowerTriMatrixViewF(
        const ConstSmallLowerTriMatrixView<T,N2,D,Si2,Sj2,C,I2>& m2) :
      mtype(m2) {}
    template <int N2, int Si2, int Sj2, IndexStyle I2>
    inline ConstLowerTriMatrixViewF(
        const SmallLowerTriMatrixView<T,N2,D,Si2,Sj2,C,I2>& m2) :
      mtype(m2) {}
    inline ~ConstLowerTriMatrixViewF() {}

  private :
    inline void operator=(const type& m2);

  }; // ConstLowerTriMatrixViewF

  template <class T, DiagType D, int Si, int Sj, bool C, IndexStyle I>
  struct Traits<LowerTriMatrixView<T,D,Si,Sj,C,I> >
  {
    typedef T value_type;

    typedef typename Traits<T>::real_type real_type;
    typedef typename Traits<T>::complex_type complex_type;
    enum { misreal = Traits<T>::isreal };
    enum { miscomplex = Traits<T>::iscomplex };

    typedef LowerTriMatrixView<T,D,Si,Sj,C,I> type;
    typedef const ConstLowerTriMatrixView<T,D,Si,Sj,C,I> calc_type;
    typedef calc_type eval_type;
    typedef LowerTriMatrix<T,D,Sj==1?RowMajor:ColMajor,I> copy_type;
    typedef InvalidType inverse_type;

    enum { mcolsize = UNKNOWN };
    enum { mrowsize = UNKNOWN };
    enum { msize = UNKNOWN };
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
    enum { mshape = munit ? UnitLowerTri : LowerTri };

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

    typedef ConstLowerTriMatrixView<T,D,mstepi,mstepj,C,I> const_view_type;
    typedef ConstLowerTriMatrixView<T,D,mstepi,mstepj,C,CStyle> const_cview_type;
    typedef ConstLowerTriMatrixView<T,D,mstepi,mstepj,C,FortranStyle> const_fview_type;
    typedef ConstLowerTriMatrixView<T,D,UNKNOWN,UNKNOWN,C> const_xview_type;
    typedef ConstLowerTriMatrixView<T,D,1,mstepj,C,I> const_cmview_type;
    typedef ConstLowerTriMatrixView<T,D,mstepi,1,C,I> const_rmview_type;
    typedef ConstLowerTriMatrixView<T,D,mstepi,mstepj,notC,I> const_conjugate_type;
    typedef ConstUpperTriMatrixView<T,D,mstepj,mstepi,C,I> const_transpose_type;
    typedef ConstUpperTriMatrixView<T,D,mstepj,mstepi,notC,I> const_adjoint_type;

    typedef ConstLowerTriMatrixView<T,NonUnitDiag,mstepi,1,C,I> const_offdiag_type;
    typedef ConstLowerTriMatrixView<T,UnitDiag,mstepi,1,C,I> const_unitdiag_type;
    typedef ConstLowerTriMatrixView<real_type,D,twoSi,twoSj,false,I> const_realview_type;
    typedef const_realview_type const_imagview_type;
    typedef ConstLowerTriMatrixView<T,D,mstepi,mstepj,false,I> const_nonconj_type;
    typedef LowerTriMatrixView<T,D,mstepi,mstepj,C,I> nonconst_type;

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

    typedef LowerTriMatrixView<T,D,mstepi,mstepj,C,I> view_type;
    typedef LowerTriMatrixView<T,D,mstepi,mstepj,C,CStyle> cview_type;
    typedef LowerTriMatrixView<T,D,mstepi,mstepj,C,FortranStyle> fview_type;
    typedef LowerTriMatrixView<T,D,UNKNOWN,UNKNOWN,C> xview_type;
    typedef LowerTriMatrixView<T,D,1,mstepj,C,I> cmview_type;
    typedef LowerTriMatrixView<T,D,mstepi,1,C,I> rmview_type;
    typedef LowerTriMatrixView<T,D,mstepi,mstepj,notC,I> conjugate_type;
    typedef UpperTriMatrixView<T,D,mstepj,mstepi,C,I> transpose_type;
    typedef UpperTriMatrixView<T,D,mstepj,mstepi,notC,I> adjoint_type;

    typedef LowerTriMatrixView<T,NonUnitDiag,mstepi,mstepj,C,I> offdiag_type;
    typedef LowerTriMatrixView<T,UnitDiag,mstepi,mstepj,C,I> unitdiag_type;
    typedef LowerTriMatrixView<real_type,D,twoSi,twoSj,false,I> realview_type;
    typedef realview_type imagview_type;
    typedef LowerTriMatrixView<T,D,mstepi,mstepj,false,I> nonconj_type;
  };

  template <class T, DiagType D, int Si, int Sj, bool C, IndexStyle I>
  class LowerTriMatrixView :
    public BaseMatrix_Tri_Mutable<LowerTriMatrixView<T,D,Si,Sj,C,I> >
  {
  public:

    typedef LowerTriMatrixView<T,D,Si,Sj,C,I> type;
    typedef BaseMatrix_Tri_Mutable<type> base_mut;
    typedef typename base_mut::reference reference;
    enum { mrowmajor = Traits<type>::mrowmajor };
    enum { mcolmajor = Traits<type>::mcolmajor };
    enum { mshape = Traits<type>::mshape };
    enum { munit = Traits<type>::munit };

    //
    // Constructors
    //

    inline LowerTriMatrixView(T* m, size_t s, int si, int sj) :
      itsm(m), itss(s), itssi(si), itssj(sj) {}

    inline LowerTriMatrixView(T* m, size_t s, int si) :
      itsm(m), itss(s), itssi(si), itssj(Sj)
    { TMVStaticAssert(Sj != UNKNOWN); }

    inline LowerTriMatrixView(T* m, size_t s) :
      itsm(m), itss(s), itssi(Si), itssj(Sj)
    { TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); }

    inline LowerTriMatrixView(const type& m2) :
      itsm(m2.itsm), itss(m2.size()), 
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int Si2, int Sj2, IndexStyle I2>
    inline LowerTriMatrixView(LowerTriMatrixView<T,D,Si2,Sj2,C,I2> m2) :
      itsm(m2.ptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int N2, int Si2, int Sj2, IndexStyle I2>
    inline LowerTriMatrixView(SmallLowerTriMatrixView<T,N2,D,Si2,Sj2,C,I2> m2) :
      itsm(m2.ptr()), itss(m2.size()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    inline ~LowerTriMatrixView() {
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
    const size_t itss;
    const StepInt<Si> itssi;
    const StepInt<Sj> itssj;

  }; // LowerTriMatrixView

  template <class T, DiagType D, int Si, int Sj, bool C>
  class LowerTriMatrixViewF :
    public LowerTriMatrixView<T,D,Si,Sj,C,FortranStyle>
  {
  public:

    typedef LowerTriMatrixViewF<T,D,Si,Sj,C> type;
    typedef LowerTriMatrixView<T,D,Si,Sj,C,FortranStyle> mtype;

    inline LowerTriMatrixViewF(T* m, size_t s, int si, int sj) :
      mtype(m,s,si,sj) {}
    inline LowerTriMatrixViewF(T* m, size_t s, int si) :
      mtype(m,s,si) {}
    inline LowerTriMatrixViewF(T* m, size_t s) : mtype(m,s) {}
    inline LowerTriMatrixViewF(const type& m2) : mtype(m2) {}
    template <int Si2, int Sj2, IndexStyle I2>
    inline LowerTriMatrixViewF(LowerTriMatrixView<T,D,Si2,Sj2,C,I2> m2) :
      mtype(m2) {}
    template <int N2, int Si2, int Sj2, IndexStyle I2>
    inline LowerTriMatrixViewF(
        SmallLowerTriMatrixView<T,N2,D,Si2,Sj2,C,I2> m2) : mtype(m2) {}
    inline ~LowerTriMatrixViewF() {}

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

  }; // LowerTriMatrixViewF



  //---------------------------------------------------------------------------

  //
  // Special Creators: 
  //   UpperTriMatrixViewOf(T* m, n, S)
  //   UpperTriMatrixViewOf(T* m, n, si, sj)
  //   UnitUpperTriMatrixViewOf(T* m, n, S)
  //   UnitUpperTriMatrixViewOf(T* m, n, si, sj)
  //   LowerTriMatrixViewOf(T* m, n, S)
  //   LowerTriMatrixViewOf(T* m, n, si, sj)
  //   UnitLowerTriMatrixViewOf(T* m, n, S)
  //   UnitLowerTriMatrixViewOf(T* m, n, si, sj)
  //

  template <class T> 
  inline UpperTriMatrixView<T,NonUnitDiag> UpperTriMatrixViewOf(
      T* m, size_t size, StorageType stor)
  {
    TMVAssert(stor == RowMajor || stor == ColMajor);
    if (stor == RowMajor)
      return UpperTriMatrixView<T,NonUnitDiag>(m,size,size,1);
    else
      return UpperTriMatrixView<T,NonUnitDiag>(m,size,1,size);
  }

  template <class T> 
  inline ConstUpperTriMatrixView<T,NonUnitDiag> UpperTriMatrixViewOf(
      const T* m, size_t size, StorageType stor)
  {
    TMVAssert(stor == RowMajor || stor == ColMajor);
    if (stor == RowMajor)
      return ConstUpperTriMatrixView<T,NonUnitDiag>(m,size,size,1);
    else
      return ConstUpperTriMatrixView<T,NonUnitDiag>(m,size,1,size);
  }

  template <class T> 
  inline UpperTriMatrixView<T,NonUnitDiag> UpperTriMatrixViewOf(
      T* m, size_t size, int stepi, int stepj)
  { return UpperTriMatrixView<T,NonUnitDiag>(m,size,stepi,stepj); }

  template <class T> 
  inline ConstUpperTriMatrixView<T,NonUnitDiag> UpperTriMatrixViewOf(
      const T* m, size_t size, int stepi, int stepj)
  { return ConstUpperTriMatrixView<T,NonUnitDiag>(m,size,stepi,stepj); }

  template <class T> 
  inline UpperTriMatrixView<T,UnitDiag> UnitUpperTriMatrixViewOf(
      T* m, size_t size, StorageType stor)
  {
    TMVAssert(stor == RowMajor || stor == ColMajor);
    if (stor == RowMajor)
      return UpperTriMatrixView<T,UnitDiag>(m,size,size,1);
    else
      return UpperTriMatrixView<T,UnitDiag>(m,size,1,size);
  }

  template <class T> 
  inline ConstUpperTriMatrixView<T,UnitDiag> UnitUpperTriMatrixViewOf(
      const T* m, size_t size, StorageType stor)
  {
    TMVAssert(stor == RowMajor || stor == ColMajor);
    if (stor == RowMajor)
      return ConstUpperTriMatrixView<T,UnitDiag>(m,size,size,1);
    else
      return ConstUpperTriMatrixView<T,UnitDiag>(m,size,1,size);
  }

  template <class T> 
  inline UpperTriMatrixView<T,UnitDiag> UnitUpperTriMatrixViewOf(
      T* m, size_t size, int stepi, int stepj)
  { return UpperTriMatrixView<T,UnitDiag>(m,size,stepi,stepj); }

  template <class T> 
  inline ConstUpperTriMatrixView<T,UnitDiag> UnitUpperTriMatrixViewOf(
      const T* m, size_t size, int stepi, int stepj)
  { return ConstUpperTriMatrixView<T,UnitDiag>(m,size,stepi,stepj); }

  template <class T> 
  inline LowerTriMatrixView<T,NonUnitDiag> LowerTriMatrixViewOf(
      T* m, size_t size, StorageType stor)
  {
    TMVAssert(stor == RowMajor || stor == ColMajor);
    if (stor == RowMajor)
      return LowerTriMatrixView<T,NonUnitDiag>(m,size,size,1);
    else
      return LowerTriMatrixView<T,NonUnitDiag>(m,size,1,size);
  }

  template <class T> 
  inline ConstLowerTriMatrixView<T,NonUnitDiag> LowerTriMatrixViewOf(
      const T* m, size_t size, StorageType stor)
  {
    TMVAssert(stor == RowMajor || stor == ColMajor);
    if (stor == RowMajor)
      return ConstLowerTriMatrixView<T,NonUnitDiag>(m,size,size,1);
    else
      return ConstLowerTriMatrixView<T,NonUnitDiag>(m,size,1,size);
  }

  template <class T> 
  inline LowerTriMatrixView<T,NonUnitDiag> LowerTriMatrixViewOf(
      T* m, size_t size, int stepi, int stepj)
  { return LowerTriMatrixView<T,NonUnitDiag>(m,size,stepi,stepj); }

  template <class T> 
  inline ConstLowerTriMatrixView<T,NonUnitDiag> LowerTriMatrixViewOf(
      const T* m, size_t size, int stepi, int stepj)
  { return ConstLowerTriMatrixView<T,NonUnitDiag>(m,size,stepi,stepj); }

  template <class T> 
  inline LowerTriMatrixView<T,UnitDiag> UnitLowerTriMatrixViewOf(
      T* m, size_t size, StorageType stor)
  {
    TMVAssert(stor == RowMajor || stor == ColMajor);
    if (stor == RowMajor)
      return LowerTriMatrixView<T,UnitDiag>(m,size,size,1);
    else
      return LowerTriMatrixView<T,UnitDiag>(m,size,1,size);
  }

  template <class T> 
  inline ConstLowerTriMatrixView<T,UnitDiag> UnitLowerTriMatrixViewOf(
      const T* m, size_t size, StorageType stor)
  {
    TMVAssert(stor == RowMajor || stor == ColMajor);
    if (stor == RowMajor)
      return ConstLowerTriMatrixView<T,UnitDiag>(m,size,size,1);
    else
      return ConstLowerTriMatrixView<T,UnitDiag>(m,size,1,size);
  }

  template <class T> 
  inline LowerTriMatrixView<T,UnitDiag> UnitLowerTriMatrixViewOf(
      T* m, size_t size, int stepi, int stepj)
  { return LowerTriMatrixView<T,UnitDiag>(m,size,stepi,stepj); }

  template <class T> 
  inline ConstLowerTriMatrixView<T,UnitDiag> UnitLowerTriMatrixViewOf(
      const T* m, size_t size, int stepi, int stepj)
  { return ConstLowerTriMatrixView<T,UnitDiag>(m,size,stepi,stepj); }



  //
  // Swap
  //

  template <class T, DiagType D, StorageType S1, IndexStyle I1, int Si2, int Sj2, bool C2, IndexStyle I2>
  inline void Swap(
      UpperTriMatrix<T,D,S1,I1>& m1, UpperTriMatrixView<T,D,Si2,Sj2,C2,I2> m2)
  { DoSwap(m1,m2); }
  template <class T, DiagType D, int Si1, int Sj1, bool C1, IndexStyle I1, StorageType S2, IndexStyle I2>
  inline void Swap(
      UpperTriMatrixView<T,D,Si1,Sj1,C1,I1> m1, UpperTriMatrix<T,D,S2,I2>& m2)
  { DoSwap(m1,m2); }
  template <class T, DiagType D, int Si1, int Sj1, bool C1, IndexStyle I1, int Si2, int Sj2, bool C2, IndexStyle I2>
  inline void Swap(
      UpperTriMatrixView<T,D,Si1,Sj1,C1,I1> m1, UpperTriMatrixView<T,D,Si2,Sj2,C2,I2> m2)
  { DoSwap(m1,m2); }

  template <class T, DiagType D, StorageType S1, IndexStyle I1, int Si2, int Sj2, bool C2, IndexStyle I2>
  inline void Swap(
      LowerTriMatrix<T,D,S1,I1>& m1, LowerTriMatrixView<T,D,Si2,Sj2,C2,I2> m2)
  { Swap(m1.Transpose(),m2.Transpose()); }
  template <class T, DiagType D, int Si1, int Sj1, bool C1, IndexStyle I1, StorageType S2, IndexStyle I2>
  inline void Swap(
      LowerTriMatrixView<T,D,Si1,Sj1,C1,I1> m1, LowerTriMatrix<T,D,S2,I2>& m2)
  { Swap(m1.Transpose(),m2.Transpose()); }
  template <class T, DiagType D, int Si1, int Sj1, bool C1, IndexStyle I1, int Si2, int Sj2, bool C2, IndexStyle I2>
  inline void Swap(
      LowerTriMatrixView<T,D,Si1,Sj1,C1,I1> m1, LowerTriMatrixView<T,D,Si2,Sj2,C2,I2> m2)
  { Swap(m1.Transpose(),m2.Transpose()); }


  //
  // TypeText 
  //

  template <class T, DiagType D, StorageType S, IndexStyle I>
  inline std::string TypeText(const UpperTriMatrix<T,D,S,I>& )
  {
    std::ostringstream s;
    s << "UpperTriMatrix<"<<TypeText(T())<<","<<Text(D);
    s << ","<<Text(S)<<","<<Text(I)<<">";
    return s.str();
  }

  template <class T, DiagType D, int Si, int Sj, bool C, IndexStyle I>
  inline std::string TypeText(const ConstUpperTriMatrixView<T,D,Si,Sj,C,I>& m)
  {
    std::ostringstream s;
    s << "ConstUpperTriMatrixView<"<<TypeText(T())<<","<<Text(D);
    s << ","<<IntTraits<Si>::text();
    if (Si == UNKNOWN) s << "("<<m.stepi()<<")";
    s << ","<<IntTraits<Sj>::text();
    if (Sj == UNKNOWN) s << "("<<m.stepj()<<")";
    s << ","<<C<<","<<Text(I)<<">";
    return s.str();
  }

  template <class T, DiagType D, int Si, int Sj, bool C, IndexStyle I>
  inline std::string TypeText(const UpperTriMatrixView<T,D,Si,Sj,C,I>& m)
  {
    std::ostringstream s;
    s << "UpperTriMatrixView<"<<TypeText(T())<<","<<Text(D);
    s << ","<<IntTraits<Si>::text();
    if (Si == UNKNOWN) s << "("<<m.stepi()<<")";
    s << ","<<IntTraits<Sj>::text();
    if (Sj == UNKNOWN) s << "("<<m.stepj()<<")";
    s << ","<<C<<","<<Text(I)<<">";
    return s.str();
  }

  template <class T, DiagType D, StorageType S, IndexStyle I>
  inline std::string TypeText(const LowerTriMatrix<T,D,S,I>& )
  {
    std::ostringstream s;
    s << "LowerTriMatrix<"<<TypeText(T())<<","<<Text(D);
    s << ","<<Text(S)<<","<<Text(I)<<">";
    return s.str();
  }

  template <class T, DiagType D, int Si, int Sj, bool C, IndexStyle I>
  inline std::string TypeText(const ConstLowerTriMatrixView<T,D,Si,Sj,C,I>& m)
  {
    std::ostringstream s;
    s << "ConstLowerTriMatrixView<"<<TypeText(T())<<","<<Text(D);
    s << ","<<IntTraits<Si>::text();
    if (Si == UNKNOWN) s << "("<<m.stepi()<<")";
    s << ","<<IntTraits<Sj>::text();
    if (Sj == UNKNOWN) s << "("<<m.stepj()<<")";
    s << ","<<C<<","<<Text(I)<<">";
    return s.str();
  }

  template <class T, DiagType D, int Si, int Sj, bool C, IndexStyle I>
  inline std::string TypeText(const LowerTriMatrixView<T,D,Si,Sj,C,I>& m)
  {
    std::ostringstream s;
    s << "LowerTriMatrixView<"<<TypeText(T())<<","<<Text(D);
    s << ","<<IntTraits<Si>::text();
    if (Si == UNKNOWN) s << "("<<m.stepi()<<")";
    s << ","<<IntTraits<Sj>::text();
    if (Sj == UNKNOWN) s << "("<<m.stepj()<<")";
    s << ","<<C<<","<<Text(I)<<">";
    return s.str();
  }


} // namespace tmv

#endif
