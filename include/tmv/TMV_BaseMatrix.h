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
// This file defines the BaseMatrix, BaseMatrix_Calc, and BaseMatrix_Mutable
// classes.  These define the functionality for matrices of all the
// different shapes, including dense rectangular, upper/lower triangular,
// diagonal, banded, symmetric, hermitian, band-symmetric and band-hermitian.
//
// BaseMatrix is the base class for all of the various types of matrices.
// It defines all of the functions that you can use on any matrix.
//
// There is a single template argument, which is the type of the most
// derived class for the object.
// So, for example, Matrix<double> inherits from BaseMatrix<Matrix<double> >.
//
// The methods that are defined for BaseMatrix are described in 
// TMV_Matrix.h and include all of the constant access methods,
// the non-modifying functions, write, and arithmetic operators.
//
//
// BaseMatrix_Calc is the base class for all matrices that have their
// data calculated in memory, possibly being the conjugate of the 
// underlying data.
//
// The methods that are defined for BaseMatrix_Calc are described in
// TMV_Matrix.h and include row, col, diag, Transpose, Conjugate, 
// and Adjoint (among others).
//
//
// BaseMatrix_Mutable is the base class for all matrices with calculated
// data (like BaseMatrix_Calc) whose values are allowed to be modified.
// 
// The methods that are defined for BaseMatrix_Mutable are described in
// TMV_Matrix.h and include the non-const versions of Transpose, Conjugate,
// etc. and some functions such as Zero(), SetAllTo(), and AddToAll()
// that work for any shape matrix.  This class is also the type of the
// argument of AssignTo(). 
//
//
#ifndef TMV_BaseMatrix_H
#define TMV_BaseMatrix_H

#include <sstream>
#include "TMV_Vec.h"
#include "TMV_Shape.h"

namespace tmv {

  // BaseMatrix is the base class for all matrix classes.
  // All non-modifying functions are defined for BaseMatrix, such as
  // Norm, Trace, Det, etc.
  template <class M>
  class BaseMatrix;

  // BaseMatrix gets a lot of the relevant information about M 
  // from the Traits class: Traits<M>
  // The types are defined as typedef statements.
  // The integer (or boolean) values are defined as enum statements.
  // See TMV_Matrix.h or TMV_SmallMatrix.h for some concrete examples.
  //
  //
  //  value_type = The type of the individual elements
  //
  //  type = shorthand for the derived type 
  //
  //  calc_type = The type of the calculated version of the matrix.
  //  i.e. where all of the values are stored in memory somewhere.
  //  This is the return type of calc()
  //  Use this type when you will access each element
  //  of a composite matrix multiple times.
  //  Also use this when you will use cptr().
  //
  //  eval_type = The type of the evaluated version of the matrix
  //  This is the return type of eval()
  //  Use this type when you will access each element
  //  of a composite matrix only once.
  //
  //  copy_type = The type of a new copy of the matrix
  //  This is the return type of copy()
  //  Use this type when you need elements in the original matrix
  //  after overwriting those elements.
  //
  //  inverse_type = The type of the inverse of the matrix
  //  This is the return type of Inverse() const
  //
  //  mcolsize = column size of matrix (aka number of rows)
  //  mrowsize = row size of matrix (aka number of columns)
  //  (Use UNKNOWN if unknown at compile time)
  //
  //  mshape = The shape of the non-zero elements of the matrix
  //
  //  mfort = does the indexing use fortran style?
  //
  //  mcalc = are the element values already calculated in memory?


  // BaseMatrix_Calc is derived from BaseMatrix, and is used
  // for matrices that have their values already calculated somewhere.
  // So composite classes inherit directly from BaseMatrix, rather
  // than BaseMatrix_Calc.
  template <class M>
  class BaseMatrix_Calc;

  // BaseMatrix_Calc adds some more requirements to the Traits<M> class:
  //
  //  mconj = is the matrix the conjugate of the underlying data?
  //  mrowmajor = is the matrix RowMajor?
  //  mcolmajor = is the matrix ColMajor?
  //  mstor = an appropriate storage class for copying the matrix
  //
  //  const_view_type = return type from View() const
  //  const_cview_type = return type from CView() const
  //  const_fview_type = return type from FView() const
  //  const_xview_type = return type from XView() const
  //  const_cmview_type = return type from CMView() const
  //  const_rmview_type = return type from RMView() const
  //
  //  const_conjugate_type = return type from Conjugate() const
  //  const_transpose_type = return type from Transpose() const
  //  const_adjoint_type = return type from Adjoint() const
  //  const_realview_type = return type from Real() const
  //  const_imagview_type = return type from Imag() const
  //  const_nonconj_type = return type from NonConj() const
  //  nonconst_type = return type from NonConst() const


  // BaseMatrix_Mutable is derived from BaseMatrix_Calc, and is used
  // for matrices that are allowed to have their data modified.
  template <class M>
  class BaseMatrix_Mutable;

  // BaseMatrix_Mutable adds:
  //
  //  reference = return type of m(i,j)
  //
  //  view_type = return type from View() 
  //  cview_type = return type from CView()
  //  fview_type = return type from FView()
  //  xview_type = return type from XView()
  //  cmview_type = return type from CMView()
  //  rmview_type = return type from RMView()
  //
  //  conjugate_type = return type from Conjugate() 
  //  transpose_type = return type from Transpose()
  //  adjoint_type = return type from Adjoint()
  //  const_nonconj_type = return type from NonConj() const

  //
  // Helper functions and values:
  //

  // These helper functions check the validity of indices according
  // to whether the matrix uses CStyle or FortranStyle indexing.
  // They also update the indices to be consistent with CStyle.
  template <bool mfort>
  inline void CheckRowIndex(int& i, int m)
  { // CStyle
    TMVAssert(i >= 0 && "row index must be in matrix");
    TMVAssert(i < m && "row index must be in matrix");
  }
  template <bool mfort>
  inline void CheckColIndex(int& j, int n)
  { // CStyle
    TMVAssert(j >= 0 && "column index must be in matrix");
    TMVAssert(j < n && "column index must be in matrix");
  }
  template <>
  inline void CheckRowIndex<true>(int& i, int m)
  { // FortranStyle
    TMVAssert(i >= 1 && "row index must be in matrix");
    TMVAssert(i <= m && "row index must be in matrix");
    --i;
  }
  template <>
  inline void CheckColIndex<true>(int& j, int n)
  { // FortranStyle
    TMVAssert(j >= 1 && "column index must be in matrix");
    TMVAssert(j <= n && "column index must be in matrix");
    --j;
  }

  // Override SameStorage for Matrix objects:
  template <class M1, class M2>
  inline bool SameStorage(const BaseMatrix<M1>& v1, const BaseMatrix<M2>& m2)
  { return false; }
  template <class V1, class M2>
  inline bool SameStorage(const BaseVector<V1>& v1, const BaseMatrix<M2>& m2)
  { return false; }
  template <class M1, class V2>
  inline bool SameStorage(const BaseMatrix<M1>& v1, const BaseVector<V2>& m2)
  { return false; }
#ifndef TMV_NO_ALIAS_CHECK
  template <class V1, class M2>
  inline bool SameStorage(
      const BaseVector_Calc<V1>& v1, const BaseMatrix_Calc<M2>& m2)
  { return v1.Real().cptr() == m2.Real().cptr(); }
  template <class M1, class V2>
  inline bool SameStorage(
      const BaseMatrix_Calc<M1>& m1, const BaseVector_Calc<V2>& v2)
  { return m1.Real().cptr() == v2.Real().cptr(); }
  template <class M1, class M2>
  inline bool SameStorage(
      const BaseMatrix_Calc<M1>& m1, const BaseMatrix_Calc<M2>& m2)
  { return m1.Real().cptr() == m2.Real().cptr(); }
#endif

  template <class M1, class M2>
  inline bool ExactSameStorage(
      const BaseMatrix<M1>& m1, const BaseMatrix<M2>& m2)
  { return false; }
  template <class V1, class M2>
  inline bool ExactSameStorage(
      const BaseVector<V1>& v1, const BaseMatrix<M2>& m2)
  { return false; }
  template <class M1, class V2>
  inline bool ExactSameStorage(
      const BaseMatrix<M1>& m1, const BaseVector<V2>& v2)
  { return false; }
  // Step checks are done for specific shapes in the various files
  // for each shaped matrix.

  // This helper class helps decide calc_type for composite classes:
  // We don't define anything, since it needs to be specialized
  // differently for each shape.
  template <class T, int shape, int cs, int rs, bool rm, bool fort>
  struct MCopyHelper;

  template <class M>
  inline typename M::value_type DoTrace(const BaseMatrix<M>& m);

  // Defined in TMV_MatrixIO.h
  template <class M>
  inline void Write(std::ostream& os, const BaseMatrix_Calc<M>& m);
  template <class M>
  inline void Write(std::ostream& os, const BaseMatrix_Calc<M>& m,
      typename M::real_type thresh) ;

  template <class M> 
  class BaseMatrix
  {
  public:

    typedef M type;
    typedef typename Traits<type>::value_type value_type;
    typedef typename Traits<type>::calc_type calc_type;
    typedef typename Traits<type>::eval_type eval_type;
    typedef typename Traits<type>::copy_type copy_type;
    typedef typename Traits<type>::inverse_type inverse_type;

    enum { mcolsize = Traits<type>::mcolsize }; 
    enum { mrowsize = Traits<type>::mrowsize }; 
    enum { mshape = Traits<type>::mshape };
    enum { mfort = Traits<type>::mfort };
    enum { mcalc = Traits<type>::mcalc };

    typedef typename Traits<value_type>::real_type real_type;
    typedef typename Traits<value_type>::complex_type complex_type;
    enum { misreal = Traits<value_type>::isreal };
    enum { miscomplex = Traits<value_type>::iscomplex };

    //
    // Constructor
    //

    inline BaseMatrix() {}
    inline BaseMatrix(const BaseMatrix<type>&) {}
    inline ~BaseMatrix() {}

  private :
    inline void operator=(const BaseMatrix<type>& m2);
  public :


    //
    // Access
    //

    inline value_type operator()(int i, int j) const 
    {
      CheckRowIndex<mfort>(i,colsize());
      CheckColIndex<mfort>(j,rowsize());
      return cref(i,j);
    }


    //
    // Functions
    //

    inline value_type Trace() const
    { return tmv::DoTrace(eval()); }

    inline value_type SumElements() const
    { return calc().SumElements(); }

    inline real_type SumAbsElements() const
    { return calc().SumAbsElements(); }

    inline real_type MaxAbsElement() const 
    { return calc().MaxAbsElement(); }

    inline real_type NormSq() const
    { return calc().NormSq(); }

    inline real_type NormSq(const real_type scale) const
    { return calc().NormSq(scale); }

    inline real_type NormF() const
    { return calc().NormF(); }

    inline real_type Norm() const
    { return NormF(); }

    inline real_type Norm1() const 
    { return calc().Norm1(); }

    inline real_type NormInf() const 
    { return calc().NormInf(); }

    template <class ret_type, class F>
    inline ret_type SumElements(const F& f) const
    { return calc().SumElements(f); }

    
    // 
    // Division Functions
    // These are overridden for Matrix<>, since it has DivHelper.
    //

    inline inverse_type Inverse() const
    { return inverse_type(mat()); }

#if 0
    inline value_type Det() const
    {
      TMVStaticAssert((Sizes<mrowsize,mcolsize>::same));
      TMVAssert(colsize() == rowsize());
      return tmv::Det(mat());
    }

    inline real_type LogDet(value_type* sign=0) const
    {
      TMVStaticAssert((Sizes<mrowsize,mcolsize>::same));
      TMVAssert(colsize() == rowsize());
      return tmv::LogDet(mat(),sign);
    }

    inline bool Singular() const
    { return Det() == value_type(0); }

    inline real_type Norm2() const
    { return calc().Norm2(); }

    inline real_type Condition() const
    { return calc().Condition(); }

    template <class M2>
    inline void Inverse(BaseMatrix_Mutable<M2>& minv) const
    { tmv::DoInverse(mat(),minv.mat()); }

    template <class M2>
    inline void InverseATA(BaseMatrix_Mutable<M2>& minv) const
    { tmv::DoInverseATA(mat(),minv.mat()); }
#endif


    // 
    // I/O
    //

    inline void Write(std::ostream& os) const
    { tmv::Write(os,calc().CView()); }
    inline void Write(std::ostream& os, real_type thresh) const
    { tmv::Write(os,calc().CView(),thresh); }


    //
    // Auxilliary routines
    //

    inline const type& mat() const 
    { return *static_cast<const type*>(this); }

    inline calc_type calc() const 
    { return static_cast<calc_type>(mat()); }

    inline eval_type eval() const 
    { return static_cast<eval_type>(mat()); }

    inline copy_type copy() const 
    { return static_cast<copy_type>(mat()); }

    inline bool IsSquare() const 
    { return Sizes<mcolsize,mrowsize>::equal || (colsize() == rowsize()); }

    inline size_t nrows() const { return mat().colsize(); }
    inline size_t ncols() const { return mat().rowsize(); }

    // Note that these last function need to be defined in a more derived
    // class than this, or an infinite loop will result when compiling.

    inline size_t colsize() const { return mat().colsize(); }
    inline size_t rowsize() const { return mat().rowsize(); }

    inline value_type cref(int i, int j) const  { return mat().cref(i,j); }

    template <class M2>
    inline void AssignTo(BaseMatrix_Mutable<M2>& m2) const
    { mat().AssignTo(m2); }

  }; // BaseMatrix

  template <class M> 
  class BaseMatrix_Calc : 
    public BaseMatrix<M>
  {
  public:

    typedef M type;
    typedef BaseMatrix<M> base;

    typedef typename base::value_type value_type;
    typedef typename base::real_type real_type;

    enum { mfort = base::mfort };
    enum { miscomplex = base::miscomplex };
    enum { mconj = Traits<type>::mconj };
    enum { mrowmajor = Traits<type>::mrowmajor }; 
    enum { mcolmajor = Traits<type>::mcolmajor }; 
    enum { mstor = Traits<type>::mstor };

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

    typedef typename Traits<type>::nonconst_type nonconst_type;



    //
    // Constructor
    //

    inline BaseMatrix_Calc() {}
    inline BaseMatrix_Calc(const BaseMatrix_Calc<M>&) {}
    inline ~BaseMatrix_Calc() {}

  private:
    void operator=(const BaseMatrix_Calc<M>&);
  public:


    // All of these functions are implemented in a derived class,
    // but these are the things that all of the differently shaped
    // matrices should be able to implement in their base class.
    // e.g. BaseMatrix_Rec, BaseMatrix_Diag, etc.

    //
    // Views
    //

    inline const_view_type View() const
    { return mat().View(); }

    inline const_cview_type CView() const
    { return mat().View(); }

    inline const_fview_type FView() const
    { return mat().View(); }

    inline const_xview_type XView() const
    { return mat().View(); }

    inline const_cmview_type CMView() const
    {
      TMVAssert(mat().iscm());
      return mat().View(); 
    }

    inline const_rmview_type RMView() const
    {
      TMVAssert(mat().isrm());
      return mat().View(); 
    }

    inline const_transpose_type Transpose() const
    { return mat().Transpose(); }

    inline const_conjugate_type Conjugate() const
    { return mat().Conjugate(); }

    inline const_adjoint_type Adjoint() const
    { return mat().Adjoint(); }

    inline const_realview_type Real() const
    { return mat().Real(); }

    inline const_imagview_type Imag() const
    { TMVStaticAssert(miscomplex); return mat().Imag(); }

    inline const_realview_type NonConj() const
    { return mat().NonConj(); }

    inline nonconst_type NonConst() const
    { return mat().NonConst(); } 



    //
    // Auxilliary routines
    //

    inline bool isconj() const { return mconj; }

    inline bool isrm() const { return mat().isrm(); }
    inline bool iscm() const { return mat().iscm(); }
    inline StorageType stor() const 
    { return isrm() ? RowMajor : iscm() ? ColMajor : NoMajor; }

    inline const type& mat() const
    { return *static_cast<const type*>(this); }

    inline size_t colsize() const { return mat().colsize(); }
    inline size_t rowsize() const { return mat().rowsize(); }

  }; // BaseMatrix_Calc

  template <class M> 
  class BaseMatrix_Mutable 
  {
  public:

    typedef M type;

    typedef typename Traits<type>::value_type value_type;
    typedef typename Traits<value_type>::real_type real_type;

    enum { mcolsize = Traits<type>::mcolsize };
    enum { mrowsize = Traits<type>::mrowsize };
    enum { mfort = Traits<type>::mfort };

    typedef typename Traits<type>::const_view_type const_view_type;
    typedef typename Traits<type>::const_cview_type const_cview_type;
    typedef typename Traits<type>::const_fview_type const_fview_type;
    typedef typename Traits<type>::const_xview_type const_xview_type;
    typedef typename Traits<type>::const_cmview_type const_cmview_type;
    typedef typename Traits<type>::const_rmview_type const_rmview_type;
    typedef typename Traits<type>::const_transpose_type const_transpose_type;
    typedef typename Traits<type>::const_conjugate_type const_conjugate_type;
    typedef typename Traits<type>::const_adjoint_type const_adjoint_type;

    typedef typename Traits<type>::view_type view_type;
    typedef typename Traits<type>::cview_type cview_type;
    typedef typename Traits<type>::fview_type fview_type;
    typedef typename Traits<type>::xview_type xview_type;
    typedef typename Traits<type>::cmview_type cmview_type;
    typedef typename Traits<type>::rmview_type rmview_type;
    typedef typename Traits<type>::transpose_type transpose_type;
    typedef typename Traits<type>::conjugate_type conjugate_type;
    typedef typename Traits<type>::adjoint_type adjoint_type;

    typedef typename Traits<type>::reference reference;


    //
    // Constructor
    //

    inline BaseMatrix_Mutable() {}
    inline BaseMatrix_Mutable(const BaseMatrix_Mutable<M>&) {}
    inline ~BaseMatrix_Mutable() {}


    //
    // Access 
    //

    inline reference operator()(int i, int j)
    {
      CheckIndex<mfort>(i,colsize());
      CheckIndex<mfort>(j,rowsize());
      return ref(i,j);
    }

    //inline value_type operator()(int i, int j) const
    //{ return base_calc::operator()(i,j); }


    //
    // Op =
    //

    inline type& operator=(const BaseMatrix_Mutable<M>& m2) 
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



    //
    // Modifying Functions
    //

    inline type& Zero() 
    { return mat().Zero(); }

    inline type& SetAllTo(value_type val) 
    { return mat().SetAllTo(val); }

    inline type& AddToAll(value_type val) 
    { return mat().AddToAll(val); }

    inline type& Clip(real_type thresh) 
    { return mat().Clip(thresh); }

    template <class F>
    inline type& ApplyToAll(const F& f)
    { return mat().ApplyToAll(f); }

    inline type& ConjugateSelf() 
    { return mat().ConjugateSelf(); }


    //
    // Views
    //

    inline view_type View() 
    { return mat().View(); }

    inline cview_type CView() 
    { return mat().CView(); }

    inline fview_type FView() 
    { return mat().FView(); }

    inline xview_type XView() 
    { return mat().XView(); }

    inline cmview_type CMView() 
    {
      TMVAssert(mat().iscm());
      return mat().CMView(); 
    }

    inline rmview_type RMView() 
    {
      TMVAssert(mat().isrm());
      return mat().RMView(); 
    }

    inline transpose_type Transpose() 
    { return mat().Transpose(); }

    inline conjugate_type Conjugate() 
    { return mat().Conjugate(); }

    inline adjoint_type Adjoint() 
    { return mat().Adjiont(); }


    //
    // I/O
    //

    inline void Read(std::istream& is)
    { mat().Read(is); }


    //
    // Arithmetic
    //

    // These operators need to be here, rather than just defining
    // non-member operator*=, etc., since the argument to a non-member
    // function would have to be BaseMatrix_Mutable& (ie. a reference).
    // But then you couldn't write something like:
    // m.Transpose() += m2;
    // since the m.Transpose() function returns a view by value, which is not
    // castable to the non-const reference argument of operator+=.
    // So we define all these here with unspecified right hand sides.
    // They just have to be valid arguments to a MultEq, AddEq, etc.
    // function defined as non-member functions for various objects.

    template <class X2>
    inline type& operator+=(const X2& x2)
    { AddEq(mat(),x2); return mat(); }

    template <class X2>
    inline type& operator-=(const X2& x2)
    { SubtractEq(mat(),x2); return mat(); }

    template <class X2>
    inline type& operator*=(const X2& x2)
    { MultEq(mat(),x2); return mat(); }

    template <class X2>
    inline type& operator/=(const X2& x2)
    { DivEq(mat(),x2); return mat(); }

    template <class X2>
    inline type& operator%=(const X2& x2)
    { RDivEq(mat(),x2); return mat(); }


    //
    // Auxilliary routines
    //

    inline const type& mat() const
    { return *static_cast<const type*>(this); }
    inline type& mat()
    { return *static_cast<type*>(this); }

    inline size_t colsize() const { return mat().colsize(); }
    inline size_t rowsize() const { return mat().rowsize(); }
    inline reference ref(int i, int j) { return mat().ref(i,j); }

  }; // BaseMatrix_Mutable


  //
  // Trace
  //

  // This one is BaseMatrix, not BaseMatrix_Calc, 
  // since it should really be called with an eval() object, not calc(), 
  // since you don't need to calculate most of the elements.
  // This is also why we need DoTrace, rather than simply Trace, since
  // we have to make sure eval() is called.

  template <class M>
  static typename M::value_type DoTrace(const BaseMatrix<M>& m)
  {
    TMVStaticAssert((Sizes<M::mrowsize,M::mcolsize>::same));
    TMVAssert(m.colsize() == m.rowsize());
    enum { size = Sizes<M::mrowsize,M::mcolsize>::size };
    const int n = size == UNKNOWN ? m.colsize() : size;
    typename M::value_type sum(0);
    for (int i=0;i<n;++i) sum += m.cref(i,i);
    return sum;
  }


  //
  // Matrix ==, != Matrix
  //

  // I don't make any effort to optimize this, since it's 
  // probably only going to be used in assert statements and 
  // the like.  If anyone really wants to test equality, the 
  // better way to do it is something like 
  // if (Norm(a-b)<=1.e-10) {...}
  //
  // Also, using crefs as I do here, means that I don't have to 
  // write different versions for each kind of matrix.
  // This will work for every possible pairing.
  template <bool rm, int cs, int rs, class M1, class M2>
  struct EqMM_Helper;

  template <int cs, int rs, class M1, class M2>
  struct EqMM_Helper<true,cs,rs,M1,M2> // rm = true
  {
    static bool eq(const M1& m1, const M2& m2)
    {
      const int M = cs == UNKNOWN ? m1.colsize() : cs;
      const int N = rs == UNKNOWN ? m1.rowsize() : rs;
      for(int i=0;i<M;++i) {
        for(int j=0;j<N;++j) {
          if (m1.cref(i,j) != m2.cref(i,j)) return false;
        }
      }
      return true;
    }
  };

  template <int cs, int rs, class M1, class M2>
  struct EqMM_Helper<false,cs,rs,M1,M2> // rm = false
  {
    static bool eq(const M1& m1, const M2& m2)
    {
      const int M = cs == UNKNOWN ? m1.colsize() : cs;
      const int N = rs == UNKNOWN ? m1.rowsize() : rs;
      for(int j=0;j<N;++j) {
        for(int i=0;i<M;++i) {
          if (m1.cref(i,j) != m2.cref(i,j)) return false;
        }
      }
      return true;
    }
  };

  template <class M1, class M2>
  inline bool CallEq(
      const BaseMatrix_Calc<M1>& m1, const BaseMatrix_Calc<M2>& m2)
  {
    TMVStaticAssert((Sizes<M1::mcolsize,M2::mcolsize>::same)); 
    TMVStaticAssert((Sizes<M1::mrowsize,M2::mrowsize>::same)); 
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    enum { cs = Sizes<M1::mcolsize,M2::mcolsize>::size };
    enum { rs = Sizes<M1::mrowsize,M2::mrowsize>::size };
    enum { rm = M2::mrowmajor || (M1::mrowmajor && !M2::mcolmajor) };
    return EqMM_Helper<rm,cs,rs,M1,M2>::eq(m1.mat(),m2.mat());
  }

  template <class M1, class M2>
  inline bool operator==(
      const BaseMatrix<M1>& m1, const BaseMatrix<M2>& m2)
  { return CallEq(m1.calc().CView(),m2.calc().CView()); }

  template <class M1, class M2>
  inline bool operator!=(
      const BaseMatrix<M1>& m1, const BaseMatrix<M2>& m2)
  { return !(m1 == m2); }



  //
  // Other Functions of Matrices
  // (These just call the method mersion.)
  //

  template <class M>
  inline typename M::value_type Trace(const BaseMatrix<M>& m)
  { return m.Trace(); }

  template <class M>
  inline typename M::real_type Norm(const BaseMatrix<M>& m)
  { return m.Norm(); }

  template <class M>
  inline typename M::real_type NormF(const BaseMatrix<M>& m)
  { return m.NormF(); }

  template <class M>
  inline typename M::real_type NormSq(const BaseMatrix<M>& m)
  { return m.NormSq(); }

  template <class M>
  inline typename M::real_type Norm1(const BaseMatrix<M>& m)
  { return m.Norm1(); }

  template <class M>
  inline typename M::real_type NormInf(const BaseMatrix<M>& m)
  { return m.NormInf(); }

  template <class M>
  inline typename M::real_type MaxAbsElement(const BaseMatrix<M>& m)
  { return m.NormInf(); }

  template <class M>
  inline typename M::value_type SumElements(const BaseMatrix<M>& m)
  { return m.SumElements(); }

  template <class M>
  inline typename M::real_type SumAbsElements(const BaseMatrix<M>& m)
  { return m.SumAbsElements(); }

  template <class M>
  inline typename M::const_conjugate_type Conjugate(const BaseMatrix_Calc<M>& m)
  { return m.Conjugate(); }

  template <class M>
  inline typename M::const_transpose_type Transpose(const BaseMatrix_Calc<M>& m)
  { return m.Transpose(); }

  template <class M>
  inline typename M::const_adjoint_type Adjoint(const BaseMatrix_Calc<M>& m)
  { return m.Adjoint(); }

  template <class M>
  inline typename M::inverse_type Inverse(const BaseMatrix<M>& m)
  { return m.Inverse(); }


  //
  // TypeText 
  //

  template <class M>
  static std::string TypeText(const BaseMatrix<M>& m)
  {
    std::ostringstream s;
    s << "BaseMatrix< "<<TypeText(m.mat())<<" >";
    return s.str();
  }

  template <class M>
  static std::string TypeText(const BaseMatrix_Calc<M>& m)
  {
    std::ostringstream s;
    s << "BaseMatrix_Calc< "<<TypeText(m.mat())<<" >";
    return s.str();
  }

  template <class M>
  static std::string TypeText(const BaseMatrix_Mutable<M>& m)
  {
    std::ostringstream s;
    s << "BaseMatrix_Mutable< "<<TypeText(m.mat())<<" >";
    return s.str();
  }

} // namespace tmv

#endif
