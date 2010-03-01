///////////////////////////////////////////////////////////////////////////////
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
// TMV_Matrix.h and include row, col, diag, transpose, conjugate, 
// and adjoint (among others).
//
//
// BaseMatrix_Mutable is the base class for all matrices with calculated
// data (like BaseMatrix_Calc) whose values are allowed to be modified.
// 
// The methods that are defined for BaseMatrix_Mutable are described in
// TMV_Matrix.h and include the non-const versions of transpose, conjugate,
// etc. and some functions such as setZero(), setAllTo(), and addToAll()
// that work for any shape matrix.  This class is also the type of the
// argument of assignTo(). 
//
//
#ifndef TMV_BaseMatrix_H
#define TMV_BaseMatrix_H

#include <sstream>
#include "TMV_BaseVector.h"
#include "TMV_Shape.h"
#include "TMV_Array.h"

namespace tmv {

    // BaseMatrix is the base class for all matrix classes.
    // All non-modifying functions are defined for BaseMatrix, such as
    // norm, trace, det, etc.
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
    //  Use this when you will access elements of a composite matrix 
    //  multiple times or need all of the elements.
    //  Also use this when you will use cptr().
    //
    //  eval_type = The type of the evaluated version of the matrix
    //  This is the return type of eval()
    //  Use this when you will only access a few elements of a 
    //  composite matrix one time each.
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
    //  const_view_type = return type from view() const
    //  const_cview_type = return type from cView() const
    //  const_fview_type = return type from fView() const
    //  const_xview_type = return type from xView() const
    //  const_cmview_type = return type from cmView() const
    //  const_rmview_type = return type from rmView() const
    //
    //  const_conjugate_type = return type from conjugate() const
    //  const_transpose_type = return type from transpose() const
    //  const_adjoint_type = return type from adjoint() const
    //  const_realpart_type = return type from realPart() const
    //  const_imagpart_type = return type from imagPart() const
    //  const_nonconj_type = return type from nonConj() const
    //  nonconst_type = return type from nonConst() const


    // BaseMatrix_Mutable is used for matrices that are allowed to have 
    // their data modified.
    // It would naturally derive from BaseMatrix_Calc, except that we want
    // to avoid having to deal with diamond inheritance patterns.
    // For example BaseMatrix_Diag derives from BaseMatrix_Calc, but not
    // from BaseMatrix_Mutable, since it represents an immutable DiagMatrix.
    // But then BaseMatrix_Diag_Mutable derives from both BaseMatrix_Diag
    // and BaseMatrix_Mutable to represent a mutable DiagMatrix.
    // So, to avoid dealing with the sticky issues of diamond inheritance,
    // we do not derive BaseMatrix_Mutable from BaseMatrix_Calc.
    template <class M>
    class BaseMatrix_Mutable;

    // BaseMatrix_Mutable adds:
    //
    //  reference = return type of m(i,j)
    //
    //  view_type = return type from view() 
    //  cview_type = return type from cView()
    //  fview_type = return type from fView()
    //  xview_type = return type from xView()
    //  cmview_type = return type from cmView()
    //  rmview_type = return type from rmView()
    //
    //  conjugate_type = return type from conjugate() 
    //  transpose_type = return type from transpose()
    //  adjoint_type = return type from adjoint()
    //  const_nonconj_type = return type from nonConj() const

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
    { return v1.realPart().cptr() == m2.realPart().cptr(); }
    template <class M1, class V2>
    inline bool SameStorage(
        const BaseMatrix_Calc<M1>& m1, const BaseVector_Calc<V2>& v2)
    { return m1.realPart().cptr() == v2.realPart().cptr(); }
    template <class M1, class M2>
    inline bool SameStorage(
        const BaseMatrix_Calc<M1>& m1, const BaseMatrix_Calc<M2>& m2)
    {
        return 
            static_cast<const void*>(m1.mat().cptr()) == 
            static_cast<const void*>(m2.mat().cptr()); 
    }
#endif

    template <class M1, class M2>
    inline bool ExactSameStorage(
        const BaseMatrix<M1>& m1, const BaseMatrix<M2>& m2)
    { return false; }
    template <class M1, class M2>
    inline bool OppositeStorage(
        const BaseMatrix<M1>& m1, const BaseMatrix<M2>& m2)
    { return false; }
    // Step checks are done for specific shapes in the various files
    // for each shaped matrix.

    // This helper class determines if there is ExactSameStorage
    // or OppositeStorage for two matrices at compile time:
    template <class M1, class M2>
    struct MStepHelper
    { 
        typedef typename M1::value_type T1;
        typedef typename M2::value_type T2;
        enum { known = (
                M1::mstepi != UNKNOWN &&
                M1::mstepj != UNKNOWN &&
                M2::mstepi != UNKNOWN &&
                M2::mstepj != UNKNOWN ) };
        enum { same = (
                Traits2<T1,T2>::sametype &&
                known &&
                M1::mstepi == int(M2::mstepi) &&
                M1::mstepj == int(M2::mstepj) ) };
        enum { opp = (
                Traits2<T1,T2>::sametype &&
                known &&
                M1::mstepi == int(M2::mstepj) &&
                M1::mstepj == int(M2::mstepi) ) };
    };

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
    inline void Write(
        std::ostream& os,
        const BaseMatrix_Calc<M>& m, typename M::real_type thresh) ;

    // Defined in TMV_QuotXM.h
    template <int ix, class T, class M>
    class QuotXM;

    // Defined in TMV_Det.h
    template <class M>
    inline typename M::value_type Det(const BaseMatrix_Calc<M>& m);
    template <class M>
    inline typename M::real_type LogDet(
        const BaseMatrix_Calc<M>& m, typename M::value_type* sign=0);
    template <class M>
    inline bool IsSingular(const BaseMatrix_Calc<M>& m);

    // Defined in TMV_InvertM.h
    template <int ix, class T, class M1, class M2>
    inline void MakeInverse(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        BaseMatrix_Mutable<M2>& minv);

    template <class M1, class M2>
    inline void MakeInverseATA(
        const BaseMatrix_Calc<M1>& m1, BaseMatrix_Mutable<M2>& mata);


    template <class M> 
    class BaseMatrix
    {
    public:
        enum { mcolsize = Traits<M>::mcolsize };
        enum { mrowsize = Traits<M>::mrowsize };
        enum { mshape = Traits<M>::mshape };
        enum { mfort = Traits<M>::mfort };
        enum { mcalc = Traits<M>::mcalc };

        typedef M type;
        typedef typename Traits<M>::value_type value_type;
        typedef typename Traits<M>::calc_type calc_type;
        typedef typename Traits<M>::eval_type eval_type;
        typedef typename Traits<M>::copy_type copy_type;
        typedef typename Traits<M>::inverse_type inverse_type;

        // Derived values:
        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;
        enum { misreal = Traits<value_type>::isreal };
        enum { miscomplex = Traits<value_type>::iscomplex };

        //
        // Constructor
        //

        inline BaseMatrix() {}
        inline BaseMatrix(const BaseMatrix<M>&) {}
        inline ~BaseMatrix() {}

    private :
        inline void operator=(const BaseMatrix<M>& m2);
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

        inline value_type trace() const
        { return tmv::DoTrace(eval()); }

        inline value_type sumElements() const
        { return calc().sumElements(); }

        inline real_type sumAbsElements() const
        { return calc().sumAbsElements(); }

        inline real_type maxAbsElement() const 
        { return calc().maxAbsElement(); }

        inline real_type normSq() const
        { return calc().normSq(); }

        inline real_type normSq(const real_type scale) const
        { return calc().normSq(scale); }

        inline real_type normF() const
        { return calc().normF(); }

        inline real_type norm() const
        { return normF(); }

        inline real_type norm1() const 
        { return calc().norm1(); }

        inline real_type normInf() const 
        { return calc().normInf(); }

        template <class ret_type, class F>
        inline ret_type sumElements(const F& f) const
        { return calc().sumElements(f); }


        // 
        // Division Functions
        //

        inline inverse_type inverse() const
        { return inverse_type(real_type(1),calc()); }

        inline value_type det() const
        {
            TMVStaticAssert((Sizes<mrowsize,mcolsize>::same));
            TMVAssert(colsize() == rowsize());
            return tmv::Det(calc());
        }

        inline real_type logDet(value_type* sign=0) const
        {
            TMVStaticAssert((Sizes<mrowsize,mcolsize>::same));
            TMVAssert(colsize() == rowsize());
            return tmv::LogDet(calc(),sign);
        }

        inline bool isSingular() const
        { 
            TMVStaticAssert((Sizes<mrowsize,mcolsize>::same));
            TMVAssert(colsize() == rowsize());
            return tmv::IsSingular(calc());
        }

        inline real_type norm2() const
        { return calc().norm2(); }

        inline real_type condition() const
        { return calc().condition(); }

        template <class M2>
        inline void makeInverse(BaseMatrix_Mutable<M2>& minv) const
        { tmv::MakeInverse(Scaling<1,real_type>(),calc(),minv); }

        template <class M2>
        inline void makeInverseATA(BaseMatrix_Mutable<M2>& mata) const
        { tmv::MakeInverseATA(calc(),mata); }


        // 
        // I/O
        //

        inline void write(std::ostream& os) const
        { tmv::Write(os,calc().cView()); }
        inline void write(std::ostream& os, real_type thresh) const
        { tmv::Write(os,calc().cView(),thresh); }


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
        inline void assignTo(BaseMatrix_Mutable<M2>& m2) const
        { mat().assignTo(m2); }

        template <class M2>
        inline void newAssignTo(BaseMatrix_Mutable<M2>& m2) const
        { mat().newAssignTo(m2); }

    }; // BaseMatrix

    template <class M> 
    class BaseMatrix_Calc : public BaseMatrix<M>
    {
    public:
        enum { mcolsize = Traits<M>::mcolsize };
        enum { mrowsize = Traits<M>::mrowsize };
        enum { mshape = Traits<M>::mshape };
        enum { mfort = Traits<M>::mfort };
        enum { mcalc = Traits<M>::mcalc };
        enum { mconj = Traits<M>::mconj };
        enum { mrowmajor = Traits<M>::mrowmajor }; 
        enum { mcolmajor = Traits<M>::mcolmajor }; 
        enum { mstor = Traits<M>::mstor };

        typedef M type;
        typedef BaseMatrix<M> base;

        typedef typename base::value_type value_type;
        typedef typename base::real_type real_type;

        typedef typename Traits<M>::const_view_type const_view_type;
        typedef typename Traits<M>::const_cview_type const_cview_type;
        typedef typename Traits<M>::const_fview_type const_fview_type;
        typedef typename Traits<M>::const_xview_type const_xview_type;
        typedef typename Traits<M>::const_cmview_type const_cmview_type;
        typedef typename Traits<M>::const_rmview_type const_rmview_type;
        typedef typename Traits<M>::const_transpose_type 
            const_transpose_type;
        typedef typename Traits<M>::const_conjugate_type 
            const_conjugate_type;
        typedef typename Traits<M>::const_adjoint_type const_adjoint_type;
        typedef typename Traits<M>::const_realpart_type const_realpart_type;
        typedef typename Traits<M>::const_imagpart_type const_imagpart_type;
        typedef typename Traits<M>::const_nonconj_type const_nonconj_type;

        typedef typename Traits<M>::nonconst_type nonconst_type;



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

        inline const_view_type view() const
        { return mat().view(); }

        inline const_cview_type cView() const
        { return mat().view(); }

        inline const_fview_type fView() const
        { return mat().view(); }

        inline const_xview_type xView() const
        { return mat().view(); }

        inline const_cmview_type cmView() const
        {
            TMVAssert(mat().iscm());
            return mat().view(); 
        }

        inline const_rmview_type rmView() const
        {
            TMVAssert(mat().isrm());
            return mat().view(); 
        }

        inline const_transpose_type transpose() const
        { return mat().transpose(); }

        inline const_conjugate_type conjugate() const
        { return mat().conjugate(); }

        inline const_adjoint_type adjoint() const
        { return mat().adjoint(); }

        inline const_realpart_type realPart() const
        { return mat().realPart(); }

        inline const_imagpart_type imagPart() const
        { TMVStaticAssert(type::miscomplex); return mat().imagPart(); }

        inline const_nonconj_type nonConj() const
        { return mat().nonConj(); }

        inline nonconst_type nonConst() const
        { return mat().nonConst(); } 



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

        enum { mcolsize = Traits<M>::mcolsize };
        enum { mrowsize = Traits<M>::mrowsize };
        enum { mshape = Traits<M>::mshape };
        enum { mfort = Traits<M>::mfort };
        enum { mcalc = Traits<M>::mcalc };

        typedef M type;

        typedef typename Traits<M>::value_type value_type;
        typedef typename Traits<M>::real_type real_type;

        typedef typename Traits<M>::const_view_type const_view_type;
        typedef typename Traits<M>::const_cview_type const_cview_type;
        typedef typename Traits<M>::const_fview_type const_fview_type;
        typedef typename Traits<M>::const_xview_type const_xview_type;
        typedef typename Traits<M>::const_cmview_type const_cmview_type;
        typedef typename Traits<M>::const_rmview_type const_rmview_type;
        typedef typename Traits<M>::const_transpose_type 
            const_transpose_type;
        typedef typename Traits<M>::const_conjugate_type 
            const_conjugate_type;
        typedef typename Traits<M>::const_adjoint_type const_adjoint_type;

        typedef typename Traits<M>::view_type view_type;
        typedef typename Traits<M>::cview_type cview_type;
        typedef typename Traits<M>::fview_type fview_type;
        typedef typename Traits<M>::xview_type xview_type;
        typedef typename Traits<M>::cmview_type cmview_type;
        typedef typename Traits<M>::rmview_type rmview_type;
        typedef typename Traits<M>::transpose_type transpose_type;
        typedef typename Traits<M>::conjugate_type conjugate_type;
        typedef typename Traits<M>::adjoint_type adjoint_type;

        typedef typename Traits<M>::reference reference;


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
            m2.assignTo(mat());
            return mat(); 
        }

        template <class M2>
        inline type& operator=(const BaseMatrix<M2>& m2) 
        {
            TMVStaticAssert((Sizes<mcolsize,M2::mcolsize>::same));
            TMVStaticAssert((Sizes<mrowsize,M2::mrowsize>::same));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            m2.assignTo(mat());
            return mat(); 
        }



        //
        // Modifying Functions
        //

        inline type& setZero() 
        { return mat().setZero(); }

        inline type& setAllTo(value_type val) 
        { return mat().setAllTo(val); }

        inline type& addToAll(value_type val) 
        { return mat().addToAll(val); }

        inline type& clip(real_type thresh) 
        { return mat().clip(thresh); }

        template <class F>
        inline type& applyToAll(const F& f)
        { return mat().applyToAll(f); }

        inline type& conjugateSelf() 
        { return mat().conjugateSelf(); }


        //
        // Views
        //

        inline view_type view() 
        { return mat().view(); }

        inline cview_type cView() 
        { return mat().cView(); }

        inline fview_type fView() 
        { return mat().fView(); }

        inline xview_type xView() 
        { return mat().xView(); }

        inline cmview_type cmView() 
        {
            TMVAssert(mat().iscm());
            return mat().cmView(); 
        }

        inline rmview_type rmView() 
        {
            TMVAssert(mat().isrm());
            return mat().rmView(); 
        }

        inline transpose_type transpose() 
        { return mat().transpose(); }

        inline conjugate_type conjugate() 
        { return mat().conjugate(); }

        inline adjoint_type adjoint() 
        { return mat().adjiont(); }


        //
        // I/O
        //

        inline void read(std::istream& is)
        { mat().read(is); }


        //
        // Arithmetic
        //

        // These operators need to be here, rather than just defining
        // non-member operator*=, etc., since the argument to a non-member
        // function would have to be BaseMatrix_Mutable& (ie. a reference).
        // But then you couldn't write something like:
        // m.transpose() += m2;
        // since the m.transpose() function returns a view by value, which is
        // not castable to the non-const reference argument of operator+=.
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
        inline value_type cref(int i, int j) const  { return mat().cref(i,j); }

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
        const int size = Sizes<M::mrowsize,M::mcolsize>::size;
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
        const int cs = Sizes<M1::mcolsize,M2::mcolsize>::size;
        const int rs = Sizes<M1::mrowsize,M2::mrowsize>::size;
        const bool rm = M2::mrowmajor || (M1::mrowmajor && !M2::mcolmajor);
        return EqMM_Helper<rm,cs,rs,M1,M2>::eq(m1.mat(),m2.mat());
    }

    template <class M1, class M2>
    inline bool operator==(
        const BaseMatrix<M1>& m1, const BaseMatrix<M2>& m2)
    { return CallEq(m1.calc().cView(),m2.calc().cView()); }

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
    { return m.trace(); }

    template <class M>
    inline typename M::real_type Norm(const BaseMatrix<M>& m)
    { return m.norm(); }

    template <class M>
    inline typename M::real_type NormF(const BaseMatrix<M>& m)
    { return m.normF(); }

    template <class M>
    inline typename M::real_type NormSq(const BaseMatrix<M>& m)
    { return m.normSq(); }

    template <class M>
    inline typename M::real_type Norm1(const BaseMatrix<M>& m)
    { return m.norm1(); }

    template <class M>
    inline typename M::real_type NormInf(const BaseMatrix<M>& m)
    { return m.normInf(); }

    template <class M>
    inline typename M::real_type MaxAbsElement(const BaseMatrix<M>& m)
    { return m.normInf(); }

    template <class M>
    inline typename M::value_type SumElements(const BaseMatrix<M>& m)
    { return m.sumElements(); }

    template <class M>
    inline typename M::real_type SumAbsElements(const BaseMatrix<M>& m)
    { return m.sumAbsElements(); }

    template <class M>
    inline typename M::const_conjugate_type Conjugate(
        const BaseMatrix_Calc<M>& m)
    { return m.conjugate(); }

    template <class M>
    inline typename M::const_transpose_type Transpose(
        const BaseMatrix_Calc<M>& m)
    { return m.transpose(); }

    template <class M>
    inline typename M::const_adjoint_type Adjoint(const BaseMatrix_Calc<M>& m)
    { return m.adjoint(); }

    template <class M>
    inline typename M::inverse_type Inverse(const BaseMatrix<M>& m)
    { return m.inverse(); }


    //
    // TMV_Text 
    //

    template <class M>
    static std::string TMV_Text(const BaseMatrix<M>& m)
    {
        std::ostringstream s;
        s << "BaseMatrix< "<<TMV_Text(m.mat())<<" >";
        return s.str();
    }

    template <class M>
    static std::string TMV_Text(const BaseMatrix_Calc<M>& m)
    {
        std::ostringstream s;
        s << "BaseMatrix_Calc< "<<TMV_Text(m.mat())<<" >";
        return s.str();
    }

    template <class M>
    static std::string TMV_Text(const BaseMatrix_Mutable<M>& m)
    {
        std::ostringstream s;
        s << "BaseMatrix_Mutable< "<<TMV_Text(m.mat())<<" >";
        return s.str();
    }

} // namespace tmv

#endif
